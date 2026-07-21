#include "lgad_tools.hh"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <regex>
#include <sstream>

using namespace Garfield;

double ElapsedS(const Clock::time_point& t0) {
  return std::chrono::duration<double>(Clock::now() - t0).count();
}

double ParseBiasFromFilename(const std::string& path) {
  const auto slash = path.find_last_of("/\\");
  const std::string base = slash == std::string::npos
      ? path : path.substr(slash + 1);
  static const std::regex re(R"(([0-9]+(?:\.[0-9]+)?)V)");
  std::smatch m;
  if (std::regex_search(base, m, re)) {
    return std::atof(m[1].str().c_str());
  }
  return std::numeric_limits<double>::quiet_NaN();
}

std::string FormatBias(double v) {
  if (std::isnan(v)) return "NA";
  std::ostringstream oss;
  if (v == std::floor(v)) {
    oss << static_cast<long long>(v);
  } else {
    oss << std::defaultfloat << std::setprecision(6) << v;
  }
  return oss.str();
}

double MedianOf(std::vector<std::size_t> v) {
  if (v.empty()) return 0.;
  const std::size_t n = v.size();
  std::nth_element(v.begin(), v.begin() + n / 2, v.end());
  const double hi = double(v[n / 2]);
  if (n % 2 == 1) return hi;
  std::nth_element(v.begin(), v.begin() + n / 2 - 1, v.begin() + n / 2);
  const double lo = double(v[n / 2 - 1]);
  return 0.5 * (lo + hi);
}

bool StripsInsideMap(const std::vector<Strip>& strips, double bx0,
                     double bx1) {
  bool ok = true;
  for (const auto& s : strips) {
    if ((s.centerUm + s.halfWidthUm) * 1.e-4 > bx1 ||
        (s.centerUm - s.halfWidthUm) * 1.e-4 < bx0) {
      std::cerr << s.label << " (" << s.centerUm << " um) lies outside the "
                << "map x = [" << bx0 * 1.e4 << ", " << bx1 * 1.e4
                << "] um -- wrong .sta?\n";
      ok = false;
    }
  }
  return ok;
}

void ScanValidity(Component& cmp, double x0, double x1, double y0,
                  double y1, int nx, int ny, const std::string& csvPath) {
  const auto tStart = Clock::now();
  std::ofstream fvalid(csvPath);
  fvalid << "x_um,y_um,status\n";
  std::size_t nInvalid = 0, nScanned = 0;
  for (int ix = 0; ix <= nx; ++ix) {
    const double x = x0 + (x1 - x0) * ix / nx;
    for (int iy = 0; iy <= ny; ++iy) {
      const double y = y0 + (y1 - y0) * iy / ny;
      double ex, ey, ez, v; int st; Medium* m = nullptr;
      cmp.ElectricField(x, y, 0., ex, ey, ez, v, m, st);
      ++nScanned;
      if (st != 0) {
        ++nInvalid;
        fvalid << x * 1.e4 << "," << y * 1.e4 << "," << st << "\n";
      }
    }
  }
  std::cout << "validity scan: " << nInvalid << "/" << nScanned
            << " points invalid; written to " << csvPath << std::endl;
  std::cout << "[timer] validity scan: " << ElapsedS(tStart) << " s"
            << std::endl;
}

/* 1D scan at x0 from y0 to y1: finds the active-silicon bounds and the
gain-layer peak, writes the profile to a text file. Prints its own
diagnostics (active silicon bounds, peak field, low-field NOTE, or the
"no valid drift medium" error if !valid).*/

FieldProfile ScanFieldProfile(Component& cmp, double x0, double y0,
                              double y1, int nScan,
                              const std::string& csvPath) {
  FieldProfile prof;
  std::ofstream fprof(csvPath);
  fprof << "y_um,V,Ex_Vcm,Ey_Vcm,Emag_Vcm\n";
  // sentinel-swap: yTop starts above yBot so "yTop > yBot" doubles as a
  // "first valid point found yet?" flag; see main.cpp's original comment.
  double yTop = 1., yBot = -1.;
  for (int i = 0; i <= nScan; ++i) {
    const double y = y0 + ((y1 - y0) * i) / nScan;
    double ex, ey, ez, v; int st; Medium* m = nullptr;
    cmp.ElectricField(x0, y, 0., ex, ey, ez, v, m, st);
    if (st != 0) continue;
    if (yTop > yBot) yTop = y;
    yBot = y;
    const double e = std::sqrt(ex * ex + ey * ey);
    if (e > prof.eMax) { prof.eMax = e; prof.yGain = y; }
    fprof << y * 1.e4 << "," << v << "," << ex << "," << ey << ","
          << e << "\n";
  }
  fprof.close();
  if (yTop > yBot) {
    std::cerr << "No valid drift medium found along the scan line --\n"
              << "check the region/material assignment above.\n";
    prof.valid = false;
    return prof;
  }
  prof.yTop = yTop;
  prof.yBot = yBot;
  prof.d = yBot - yTop;
  prof.valid = true;
  std::cout << "active silicon: y = [" << yTop * 1.e4 << ", " << yBot * 1.e4
            << "] um (d = " << prof.d * 1.e4 << " um)" << std::endl;
  std::cout << "Peak field " << prof.eMax << " V/cm at y = "
            << prof.yGain * 1.e4 << " um  (gain layer)" << std::endl;
  if (prof.eMax < 2.5e5) {
    std::cout << "NOTE: peak field < 250 kV/cm -- expect gain near 1 "
              << "at this bias." << std::endl;
  }
  return prof;
}

/* E-field over an nx*ny grid spanning [x0,x1]x[y0,y1], to a text file:
x_um,y_um,Ex_Vcm,Ey_Vcm,Emag_Vcm,V. 

Does not print -- callers want different summary lines (window vs full-device); use the returned result to print. */
FieldDumpResult DumpElectricField(Component& cmp, double x0, double x1,
                                  double y0, double y1, int nx, int ny,
                                  const std::string& csvPath, double eps) {
  std::ofstream fdump(csvPath);
  fdump << "x_um,y_um,Ex_Vcm,Ey_Vcm,Emag_Vcm,V\n";
  FieldDumpResult r;
  for (int ix = 0; ix <= nx; ++ix) {
    const double x = x0 + (x1 - x0) * ix / nx + eps;
    for (int iy = 0; iy <= ny; ++iy) {
      const double y = y0 + (y1 - y0) * iy / ny;
      double ex, ey, ez, v; int st; Medium* m = nullptr;
      cmp.ElectricField(x, y, 0., ex, ey, ez, v, m, st);
      if (st != 0) { ++r.nBad; continue; }
      const double e = std::sqrt(ex * ex + ey * ey);
      if (e > r.wMax) { r.wMax = e; r.wMaxX = x; r.wMaxY = y; }
      fdump << x * 1.e4 << "," << y * 1.e4 << "," << ex << "," << ey
            << "," << e << "," << v << "\n";
      ++r.nRows;
    }
  }
  return r;
}

/* weighting field
- per-strip Ew/wpot startup check (probed just inside the gap) 
-plus an overlap sanity check under strip0. Pure diagnostic printing. */
void PrintWeightingSanity(ComponentAnalyticField& wcmp,
                          const std::vector<Strip>& strips, double yTop,
                          double yBot) {
  for (const auto& s : strips) {
    const double xc = s.centerUm * 1.e-4, yMid = 0.5 * (yTop + yBot);
    double wx = 0., wy = 0., wz = 0.;
    wcmp.WeightingField(xc, yMid, 0., wx, wy, wz, s.label);
    std::cout << s.label << " (x=" << s.centerUm << " um): Ew at centre = ("
              << wx << ", " << wy << ", " << wz << "), wpot: top="
              << wcmp.WeightingPotential(xc, yTop + 0.01e-4, 0., s.label)
              << " mid=" << wcmp.WeightingPotential(xc, yMid, 0., s.label)
              << " back="
              << wcmp.WeightingPotential(xc, yBot - 0.01e-4, 0., s.label)
              << std::endl;
  }
  if (strips.empty()) return;
  // sum every strip's potential AT strip[0]'s own location -- generalises
  // the original's hardcoded strip0/1/2 sum to however many strips there are
  double sumAtFirst = 0.;
  for (const auto& s : strips) {
    sumAtFirst += wcmp.WeightingPotential(strips[0].centerUm * 1.e-4,
                                          yTop + 0.01e-4, 0., s.label);
  }
  std::cout << "sum of weighting potentials at yTop under "
            << strips[0].label << "'s centre (sanity: strips shouldn't "
            << "overlap) = " << sumAtFirst << "  [expect ~1]" << std::endl;
}

/* Avalanche Stepping
- step function (fine band around gain layer, coarse elsewhere)
- size cap
- time window
- stepCm captured by ref so RunConvergenceScan can tune it
- no `static` counters in the step lambda -- called twice (aval,
  avalLadder), so statics would be shared across both closures */
void ConfigureAvalanche(AvalancheMC& av, double& stepCm, double bulkStepCm,
                        double yFineLoCm, double yFineHiCm,
                        double timeWindowNs, std::size_t sizeCap,
                        bool enableSignal, bool enableMultithreading,
                        const std::string& tag,
                        std::atomic<long long>& nCalls,
                        std::atomic<long long>& nFine,
                        std::atomic<long long>& nCoarse,
                        std::atomic<long long>& nPrinted) {
  av.EnableMultithreading(enableMultithreading);
  av.EnableSignalCalculation(enableSignal);
  av.SetTimeWindow(0., timeWindowNs);
  av.EnableAvalancheSizeLimit(sizeCap);
  av.SetStepDistanceFunction(
      [&stepCm, bulkStepCm, yFineLoCm, yFineHiCm, tag, &nCalls, &nFine,
       &nCoarse, &nPrinted](double x, double y, double z) {
        ++nCalls;
        if (nPrinted.fetch_add(1) < 5) {
          std::cout << "[stepfn " << tag << "] raw args: (" << x * 1.e4
                    << ", " << y * 1.e4 << ", " << z * 1.e4 << ") um"
                    << std::endl;
        }
        if (y > yFineLoCm && y < yFineHiCm) { ++nFine; return stepCm; }
        ++nCoarse;
        return bulkStepCm;
      });
}

// per-strip weighting potential over a grid, to a text file:
// x_um,y_um,<label>_phi,...
void DumpWeightingField(ComponentAnalyticField& wcmp,
                        const std::vector<Strip>& strips, double bx0,
                        double bx1, double yTop, double yBot,
                        const std::string& outPath, int nx, int ny) {
  const auto tStart = Clock::now();
  std::ofstream fw(outPath);
  fw << "x_um,y_um";
  for (const auto& s : strips) fw << "," << s.label << "_phi";
  fw << "\n";
  for (int ix = 0; ix <= nx; ++ix) {
    const double x = bx0 + (bx1 - bx0) * ix / nx;
    for (int iy = 0; iy <= ny; ++iy) {
      const double y = yTop + (yBot - yTop) * iy / ny;
      fw << x * 1.e4 << "," << y * 1.e4;
      for (const auto& s : strips) {
        fw << "," << wcmp.WeightingPotential(x, y, 0., s.label);
      }
      fw << "\n";
    }
  }
  std::cout << "wrote " << outPath << " (" << (nx + 1) * (ny + 1)
            << " points)" << std::endl;
  std::cout << "[timer] weighting field dump: " << ElapsedS(tStart)
            << " s" << std::endl;
}

/* convergence scan
G_e/G_eh vs step size at a fixed injection point; avalLadder must
already be configured via ConfigureAvalanche. */
void RunConvergenceScan(AvalancheMC& avalLadder, double& stepCm, double x0,
                        double yInj, std::size_t sizeCap,
                        const std::string& outDir) {
  const double distanceStepsCm[] = {1.e-5, 5.e-6, 2.e-6, 1.e-6, 5.e-7, 1.e-7};
  std::vector<std::size_t> sizes;  // pooled only from the finest step
  std::size_t nCapped = 0;
  bool ehDivergent = false;
  std::cout << "# step[nm]   G_e (mean+-sem, N)      G_eh (mean+-sem, N)\n";
  const auto tLadderStart = Clock::now();
  for (std::size_t iStep = 0; iStep < std::size(distanceStepsCm); ++iStep) {
    const double thisStepCm = distanceStepsCm[iStep];
    const bool isFinest = (iStep + 1 == std::size(distanceStepsCm));
    stepCm = thisStepCm;  // step function captured stepCm by reference
    double res[2][2] = {{0., 0.}, {0., 0.}};  // [mode][sum, sum2]
    int nDone[2] = {0, 0};
    int nCapRung[2] = {0, 0};
    const int nWant[2] = {200, 300};  // high-stat pass
    for (int mode = 0; mode < 2; ++mode) {
      if (mode == 1 && ehDivergent) {
        std::cout << "step = " << thisStepCm * 1.e7
                  << " nm   G_eh skipped (divergent at coarser step)"
                  << std::endl;
        continue;
      }
      const auto tModeStart = Clock::now();
      for (int i = 0; i < nWant[mode]; ++i) {
        if (mode == 1 && i % 25 == 0) {
          std::cout << "  e+h injection " << i << std::endl;
        }
        const double xi = x0 + (i % 5 - 2) * 2.e-4;
        if (mode == 0) {
          avalLadder.AvalancheElectron(xi, yInj, 0., 0.);
        } else {
          avalLadder.AvalancheElectronHole(xi, yInj, 0., 0.);
        }
        std::size_t ne = 0, ni = 0;
        avalLadder.GetAvalancheSize(ne, ni);
        res[mode][0] += ne;
        res[mode][1] += double(ne) * double(ne);
        ++nDone[mode];
        if (ne >= sizeCap) ++nCapRung[mode];
        if (mode == 1 && isFinest) {
          sizes.push_back(ne);
          if (ne >= sizeCap) ++nCapped;
        }
      }
      const int n = nDone[mode];
      const double mean = res[mode][0] / n;
      const double var = res[mode][1] / n - mean * mean;
      const double sem = std::sqrt(var > 0. ? var / n : 0.);
      std::cout << "step = " << thisStepCm * 1.e7 << " nm   "
                << (mode == 0 ? "G_e  = " : "G_eh = ")
                << mean << " +- " << sem << " (N=" << n << ", capped="
                << nCapRung[mode] << ")"
                << "  [timer] " << ElapsedS(tModeStart) << " s"
                << std::endl;
      if (mode == 1 && nCapRung[1] * 10 > nWant[1]) {
        ehDivergent = true;
        std::cout << "NOTE: >10% of e+h avalanches hit the size cap ("
                  << sizeCap << ") -- hole-feedback divergence (f >= 1) "
                  << "at this bias for this ionisation model. G_eh is not "
                  << "a defined quantity here; skipping remaining e+h "
                  << "rungs. Compare models via G_e at this bias, or "
                  << "re-solve the Silvaco deck at lower bias."
                  << std::endl;
      }
    }
  }
  std::cout << "[timer] full ladder: " << ElapsedS(tLadderStart)
            << " s" << std::endl;

  if (!sizes.empty()) {
    std::size_t mx = 0, nBig = 0;
    double mean = 0., mean2 = 0.;
    for (const auto s : sizes) { mean += s; mean2 += double(s) * double(s); }
    mean /= sizes.size();
    mean2 /= sizes.size();
    for (const auto s : sizes) {
      if (s > mx) mx = s;
      if (s > 5. * mean) ++nBig;
    }
    const double excessNoiseF = mean > 0. ? mean2 / (mean * mean) : 0.;
    std::cout << "e+h avalanche tail (finest step, " << sizes.size()
              << " events): max = " << mx << ", "
              << nBig << "/" << sizes.size() << " above 5x mean, "
              << nCapped << " capped at " << sizeCap
              << " (feedback-divergent)" << std::endl;
    std::cout << "excess noise factor F = <G^2>/<G>^2 = " << excessNoiseF
              << "  (McIntyre low-k limit -> F -> 1; large F/heavy tail "
              << "suggests operation near the breakdown knee)"
              << std::endl;
    std::ofstream fs(outDir + "eh_sizes.txt");
    for (const auto s : sizes) fs << s << "\n";
    std::cout << "avalanche sizes (finest step only) written to "
              << "eh_sizes.txt" << std::endl;
  }
}

/* impact-ionisation model comparison
- G_e/G_eh for vodm/okuto/massey/grant at one step size and injection point
- writes eh/ge_sizes*.txt and appends to results.csv
- Caller must reset si's model afterward. */
void RunModelComparison(AvalancheMC& avalLadder, MediumSilicon& si,
                        double x0, double yInj, double stepCm,
                        std::size_t ladderCap, const std::string& outDir,
                        const std::string& biasLabel, double eMax,
                        double gMax) {
  const char* cmpModels[] = {"vodm", "okuto", "massey", "grant"};
  const int nCmpModels = 4;
  struct ModelResult {
    double ge = 0., geSem = 0., geMed = 0., geh = 0., gehSem = 0., F = 0.,
           med = 0.;
    int nCapEh = 0, nEh = 0;
    bool divergent = false;
  } cmpRes[4];

  const auto tCmpStart = Clock::now();
  std::cout << "# model comparison: fine step = " << stepCm * 1.e7
            << " nm, cap = " << ladderCap << std::endl;
  for (int im = 0; im < nCmpModels; ++im) {
    const std::string m = cmpModels[im];
    if (m == "massey") {
      si.SetImpactIonisationModelMassey();
    } else if (m == "okuto") {
      si.SetImpactIonisationModelOkutoCrowell();
    } else if (m == "grant") {
      si.SetImpactIonisationModelGrant();
    } else {
      si.SetImpactIonisationModelVanOverstraetenDeMan();
    }
    ModelResult& r = cmpRes[im];
    const auto tGeStart = Clock::now();
    double sum = 0., sum2 = 0.;
    const int nWantE = 200;
    std::vector<std::size_t> geSizes;
    for (int i = 0; i < nWantE; ++i) {
      const double xi = x0 + (i % 5 - 2) * 2.e-4;
      avalLadder.AvalancheElectron(xi, yInj, 0., 0.);
      std::size_t ne = 0, ni = 0;
      avalLadder.GetAvalancheSize(ne, ni);
      sum += ne;
      sum2 += double(ne) * double(ne);
      geSizes.push_back(ne);
    }
    r.ge = sum / nWantE;
    {
      const double var = sum2 / nWantE - r.ge * r.ge;
      r.geSem = std::sqrt(var > 0. ? var / nWantE : 0.);
    }
    r.geMed = MedianOf(geSizes);
    std::cout << m << "   G_e  = " << r.ge << " +- " << r.geSem
              << " (N=" << nWantE << ", median=" << r.geMed << ")  [timer] "
              << ElapsedS(tGeStart) << " s" << std::endl;

    const auto tEhStart = Clock::now();
    sum = 0.;
    sum2 = 0.;
    std::vector<std::size_t> mSizes;
    const int nWantEh = 300;
    for (int i = 0; i < nWantEh; ++i) {
      if (i % 50 == 0) std::cout << "  e+h injection " << i << std::endl;
      const double xi = x0 + (i % 5 - 2) * 2.e-4;
      avalLadder.AvalancheElectronHole(xi, yInj, 0., 0.);
      std::size_t ne = 0, ni = 0;
      avalLadder.GetAvalancheSize(ne, ni);
      sum += ne;
      sum2 += double(ne) * double(ne);
      mSizes.push_back(ne);
      if (ne >= ladderCap) ++r.nCapEh;
      const bool hardTrip = r.nCapEh >= 6 && r.nCapEh * 5 > i + 1;
      const bool softTrip = i == 49 && r.nCapEh * 10 > 50;
      if (hardTrip || softTrip) {
        r.divergent = true;
        std::cout << m << ": " << r.nCapEh << "/" << (i + 1)
                  << " e+h avalanches hit the cap -- hole-feedback "
                  << "divergence (f >= 1) at this bias; G_eh undefined, "
                  << "aborting this model's e+h pass." << std::endl;
        break;
      }
    }
    r.nEh = int(mSizes.size());
    r.geh = sum / r.nEh;
    {
      const double var = sum2 / r.nEh - r.geh * r.geh;
      r.gehSem = std::sqrt(var > 0. ? var / r.nEh : 0.);
    }
    {
      double m1 = 0., m2 = 0.;
      for (const auto s : mSizes) { m1 += s; m2 += double(s) * double(s); }
      m1 /= r.nEh;
      m2 /= r.nEh;
      r.F = m1 > 0. ? m2 / (m1 * m1) : 0.;
      r.med = MedianOf(mSizes);
    }
    if (!r.divergent && r.nCapEh * 10 > r.nEh) r.divergent = true;
    std::cout << m << "   G_eh = " << r.geh << " +- " << r.gehSem
              << " (N=" << r.nEh << ", median=" << r.med << ", capped="
              << r.nCapEh << (r.divergent ? ", DIVERGENT" : "")
              << ")   F = " << r.F << "  [timer] " << ElapsedS(tEhStart)
              << " s" << std::endl;
    std::ofstream fs(outDir + (biasLabel == "NA" ? ("eh_sizes_" + m + ".txt")
                                        : ("eh_sizes_" + m + "_" + biasLabel
                                           + "V.txt")));
    for (const auto s : mSizes) fs << s << "\n";
    std::ofstream fge(outDir + (biasLabel == "NA" ? ("ge_sizes_" + m + ".txt")
                                         : ("ge_sizes_" + m + "_" + biasLabel
                                            + "V.txt")));
    for (const auto s : geSizes) fge << s << "\n";
  }
  std::cout << "[timer] model comparison: " << ElapsedS(tCmpStart)
            << " s" << std::endl;
  std::cout << "\n# model   G_e (median)         G_eh             median   F"
            << "      status" << std::endl;
  for (int im = 0; im < nCmpModels; ++im) {
    const ModelResult& r = cmpRes[im];
    std::cout << cmpModels[im] << "   " << r.ge << " +- " << r.geSem
              << " (" << r.geMed << ")"
              << "   " << r.geh << " +- " << r.gehSem << "   " << r.med
              << "   " << r.F << "   "
              << (r.divergent ? "DIVERGENT (G_eh unreliable)" : "ok")
              << std::endl;
  }

  {
    std::ifstream ftest(outDir + "results.csv");
    const bool exists = ftest.good();
    ftest.close();
    std::ofstream fres(outDir + "results.csv", std::ios::app);
    if (!exists) {
      fres << "bias,model,Epeak_line_Vcm,Epeak_global_Vcm,Ge,GeSem,Geh,"
           << "GehSem,median,F,nCapEh,N,divergent,GeMedian\n";
    }
    for (int im = 0; im < nCmpModels; ++im) {
      const ModelResult& r = cmpRes[im];
      fres << biasLabel << "," << cmpModels[im] << "," << eMax << ","
           << gMax << "," << r.ge << "," << r.geSem << "," << r.geh << ","
           << r.gehSem << "," << r.med << "," << r.F << "," << r.nCapEh
           << "," << r.nEh << "," << (r.divergent ? 1 : 0) << ","
           << r.geMed << "\n";
    }
    std::cout << "appended " << nCmpModels << " rows to results.csv "
              << "(bias=" << biasLabel << ")" << std::endl;
  }
}