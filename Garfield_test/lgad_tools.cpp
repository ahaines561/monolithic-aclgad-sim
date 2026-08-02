#include "lgad_tools.hh"

#include "Garfield/Sensor.hh"
#include "Garfield/Random.hh"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <regex>
#include <stdexcept>
#include <sstream>

#if defined(_OPENMP)
#include <omp.h>
#endif

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

bool SetImpactIonisationModel(MediumSilicon& si,
                              const std::string& model) {
  if (model == "okuto") {
    si.SetImpactIonisationModelOkutoCrowell();
  } else if (model == "massey") {
    si.SetImpactIonisationModelMassey();
  } else if (model == "grant") {
    si.SetImpactIonisationModelGrant();
  } else if (model == "vodm") {
    si.SetImpactIonisationModelVanOverstraetenDeMan();
  } else {
    std::cerr << "Unknown impact-ionisation model '" << model
              << "'. Use okuto, massey, grant, or vodm.\n";
    return false;
  }
  return true;
}

void DumpTownsendCoefficients(MediumSilicon& si,
                              const std::string& csvPath,
                              const double eMinVcm,
                              const double eMaxVcm,
                              const double eStepVcm) {
  if (eStepVcm <= 0. || eMaxVcm < eMinVcm) {
    std::cerr << "DumpTownsendCoefficients: invalid field range.\n";
    return;
  }

  std::ofstream out(csvPath);
  if (!out) {
    std::cerr << "DumpTownsendCoefficients: could not open "
              << csvPath << "\n";
    return;
  }

  out << "E_Vcm,alpha_e_cm-1,alpha_h_cm-1,ok_e,ok_h\n";
  for (double e = eMinVcm; e <= eMaxVcm + 0.5 * eStepVcm;
       e += eStepVcm) {
    double alphaE = 0.;
    double alphaH = 0.;
    const bool okE =
        si.ElectronTownsend(0., e, 0., 0., 0., 0., alphaE);
    const bool okH =
        si.HoleTownsend(0., e, 0., 0., 0., 0., alphaH);
    out << e << "," << alphaE << "," << alphaH << ","
        << (okE ? 1 : 0) << "," << (okH ? 1 : 0) << "\n";
  }
  std::cout << "wrote Garfield Townsend table to " << csvPath << std::endl;
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

bool StripsInsideMap(const std::vector<ReadoutStrip>& strips, double bx0,
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
                  double y1, int nx, int ny, const std::string& csvPath,
                  double eps) {
  const auto tStart = Clock::now();
  std::ofstream fvalid(csvPath);
  fvalid << "x_um,y_um,status\n";
  std::size_t nInvalid = 0, nScanned = 0;
  // nudge off exact mesh coordinates: a point landing exactly on a node
  // or element edge belongs to neither adjoining element under strict
  // inequality tests, so lookup fails and the point reads as invalid even
  // where the region is perfectly good. DumpElectricField carries the
  // same eps for the same reason.
  for (int ix = 0; ix <= nx; ++ix) {
    const double x = x0 + (x1 - x0) * ix / nx + eps;
    for (int iy = 0; iy <= ny; ++iy) {
      const double y = y0 + (y1 - y0) * iy / ny + eps;
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
                              const std::string& csvPath,
                              double eMinVcm) {
  FieldProfile prof;
  std::ofstream fprof(csvPath);
  fprof << "y_um,V,Ex_Vcm,Ey_Vcm,Emag_Vcm\n";
  // sentinel-swap: yTop starts above yBot so "yTop > yBot" doubles as a
  // "first valid point found yet?" flag; see main.cpp's original comment.
  double yTop = 1., yBot = -1.;
  std::size_t nDead = 0;
  for (int i = 0; i <= nScan; ++i) {
    const double y = y0 + ((y1 - y0) * i) / nScan;
    double ex, ey, ez, v; int st; Medium* m = nullptr;
    cmp.ElectricField(x0, y, 0., ex, ey, ez, v, m, st);
    if (st != 0) continue;
    const double e = std::sqrt(ex * ex + ey * ey);
    // still dumped for diagnostics, but excluded from the active bounds:
    // carriers in a valid-but-fieldless region (undepleted top layer)
    // have no drift velocity, cannot leave, and burn the whole time
    // window in tiny steps -- the dominant runtime cost when included.
    if (e < eMinVcm) { ++nDead; }
    else {
      if (yTop > yBot) yTop = y;
      yBot = y;
      if (e > prof.eMax) { prof.eMax = e; prof.yGain = y; }
    }
    fprof << y * 1.e4 << "," << v << "," << ex << "," << ey << ","
          << e << "\n";
  }
  fprof.close();
  if (nDead > 0) {
    std::cout << "field profile: " << nDead << " valid point(s) below "
              << eMinVcm << " V/cm excluded from the active region"
              << std::endl;
  }
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
void PrintWeightingSanity(Component& wcmp,
                          const std::vector<ReadoutStrip>& strips, double yTop,
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
  // AvalancheMC has no DisableMultithreading(). The only entry point is
  // EnableMultithreading(unsigned nThreads = 0), and nThreads < 1 is
  // interpreted as "use hardware_concurrency()" (AvalancheMC.cc:1994).
  // Passing a bool therefore requested MAXIMUM threading, not none.
  // m_nThreads defaults to 1 (AvalancheMC.hh:347), so 1 == serial.
  av.EnableMultithreading(enableMultithreading ? 0u : 1u);
  // Pin the signal path. m_useWeightingPotential defaults to true
  // (AvalancheMC.hh:317), but the choice must match which weighting store is
  // actually populated: the potential path calls WeightingPotential() and
  // DelayedWeightingPotentials(); the field path calls WeightingField() and
  // DelayedWeightingField(). A mismatch is either an out-of-bounds read or a
  // silently-zero delayed signal, so do not rely on the default.
  av.UseWeightingPotential(true);
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

namespace {
// Production avalanche configuration. Unlike ConfigureAvalanche, this
// callback performs no atomic diagnostics on every drift step. The
// instrumented version remains available for convergence/feedback scans.
void ConfigureAvalancheFast(AvalancheMC& av, const double fineStepCm,
                            const double bulkStepCm,
                            const double yFineLoCm,
                            const double yFineHiCm,
                            const double timeWindowNs,
                            const std::size_t sizeCap,
                            const bool enableSignal) {
  // See ConfigureAvalanche: DisableMultithreading() does not exist; 1 thread
  // is serial, 0 would request hardware_concurrency(). Outer OpenMP
  // parallelism over primaries must not nest with Garfield's internal threads.
  av.EnableMultithreading(1u);
  av.UseWeightingPotential(true);
  av.EnableSignalCalculation(enableSignal);
  av.SetTimeWindow(0., timeWindowNs);
  av.EnableAvalancheSizeLimit(sizeCap);
  av.SetStepDistanceFunction(
      [fineStepCm, bulkStepCm, yFineLoCm, yFineHiCm]
      (double, const double y, double) noexcept {
        return y > yFineLoCm && y < yFineHiCm
            ? fineStepCm : bulkStepCm;
      });
}
}  // namespace

// per-strip weighting potential over a grid, to a text file:
// x_um,y_um,<label>_phi,...
void DumpWeightingField(Component& wcmp,
                        const std::vector<ReadoutStrip>& strips, double bx0,
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

DelayedWeightingValidation ValidateDelayedWeightingData(
    Component& weightingCmp, const std::vector<ReadoutStrip>& strips,
    const double probeYCm) {
  DelayedWeightingValidation out;
  if (strips.empty()) {
    out.message = "no readout electrodes were configured";
    return out;
  }

  for (std::size_t k = 0; k < strips.size(); ++k) {
    const auto& label = strips[k].label;
    const auto& times = weightingCmp.DelayedSignalTimes(label);
    if (times.empty()) {
      out.message = "electrode '" + label +
          "' has no delayed-weighting time slices";
      return out;
    }
    if (!std::is_sorted(times.begin(), times.end()) ||
        std::adjacent_find(times.begin(), times.end(),
                           std::greater_equal<double>()) != times.end()) {
      out.message = "electrode '" + label +
          "' has a delayed-weighting time grid that is not strictly increasing";
      return out;
    }
    if (times.front() < 0.) {
      out.message = "electrode '" + label +
          "' has a negative delayed-weighting time";
      return out;
    }
    if (k == 0) {
      out.timesNs = times;
    } else if (times.size() != out.timesNs.size() ||
               !std::equal(times.begin(), times.end(), out.timesNs.begin(),
                           [](const double a, const double b) {
                             const double scale =
                                 std::max({1., std::abs(a), std::abs(b)});
                             return std::abs(a - b) <= 1.e-10 * scale;
                           })) {
      out.message = "electrode '" + label +
          "' uses a different delayed-weighting time grid";
      out.timesNs.clear();
      return out;
    }
  }

  // ---------------------------------------------------------------------
  // Live probe of the store the signal path will actually read.
  //
  // The checks above only inspect DelayedSignalTimes(), which is NOT
  // sufficient: ComponentTcadBase::DelayedSignalTimes falls back from the
  // potential store (m_dwtp) to the field store (m_dwtf), so it returns a
  // non-empty grid in BOTH mismatch directions. We use the weighting-potential
  // signal path (AvalancheMC::UseWeightingPotential(true)), which consumes
  // WeightingPotential() and DelayedWeightingPotentials(). If only dynamic
  // *fields* were loaded, the delayed signal is identically zero and the only
  // warning goes to stderr -- which silenceGarfieldStderr sends to /dev/null.
  // So probe the potentials directly and fail loudly here instead.
  if (probeYCm > 0.) {
    for (const auto& s : strips) {
      const double xCm = s.centerUm * 1.e-4;
      const double wpPrompt =
          weightingCmp.WeightingPotential(xCm, probeYCm, 0., s.label);
      if (!std::isfinite(wpPrompt)) {
        out.message = "electrode '" + s.label +
            "' returned a non-finite prompt weighting potential";
        out.timesNs.clear();
        return out;
      }
      // Pick a probe time strictly inside the grid and away from t = 0.
      // At t = 0 the delayed remainder is zero BY CONSTRUCTION (Garfield
      // subtracts the prompt component), so probing there would report a
      // false failure. ComponentTcad2d never stores t = 0 in the delayed
      // grid, but ComponentGrid and ComponentUser can.
      double tMid = -1.;
      for (std::size_t i = out.timesNs.size(); i-- > 0;) {
        if (out.timesNs[i] > 0.) tMid = out.timesNs[i];
        if (i <= out.timesNs.size() / 2 && tMid > 0.) break;
      }
      if (tMid <= 0.) {
        out.message = "electrode '" + s.label +
            "' has no delayed-weighting time greater than zero";
        out.timesNs.clear();
        return out;
      }
      double dwpMax = 0.;
      for (const double y : {probeYCm, 0.5 * probeYCm, 2. * probeYCm}) {
        const double dwp =
            weightingCmp.DelayedWeightingPotential(xCm, y, 0., tMid, s.label);
        if (!std::isfinite(dwp)) {
          out.message = "electrode '" + s.label +
              "' returned a non-finite delayed weighting potential at t = " +
              std::to_string(tMid) + " ns";
          out.timesNs.clear();
          return out;
        }
        dwpMax = std::max(dwpMax, std::abs(dwp));
      }
      if (dwpMax <= 0.) {
        out.message = "electrode '" + s.label +
            "' has a delayed-weighting time grid but its delayed weighting "
            "POTENTIAL is identically zero at t = " + std::to_string(tMid) +
            " ns. The potential store is empty -- dynamic weighting FIELDS "
            "were probably loaded instead of potentials, which the "
            "UseWeightingPotential(true) signal path cannot read.";
        out.timesNs.clear();
        return out;
      }
      out.probedDelayedMax = std::max(out.probedDelayedMax, dwpMax);
      static bool printedProbe = false;
      if (printedProbe) continue;
      std::cout << "[weighting probe] " << s.label
                << ": prompt wp = " << wpPrompt << ", max |delayed wp| at t = "
                << tMid << " ns = " << dwpMax << std::endl;
      if (&s == &strips.back()) printedProbe = true;
    }
  } else {
    std::cout << "[weighting probe] SKIPPED (no probe depth given); "
              << "time-grid checks only" << std::endl;
  }

  out.available = true;
  out.message = "dynamic weighting data available and non-null for all "
                "electrodes";
  return out;
}

/* parallel avalanche pass (OpenMP)
- one Sensor + one AvalancheMC per thread; only the read-only field
  components are shared. Mirrors the pattern proven in 1000MIPs.cpp.
- per-thread signals are summed at the end (induced signals are additive)
- without -fopenmp the pragmas are ignored and this runs serially */
unsigned long RunAvalanchePass(
    Component& driftCmp, Component& weightingCmp,
    const std::vector<Component*>& perStripWeighting,
    const std::vector<ReadoutStrip>& strips,
    const std::vector<std::array<double, 4>>& primaries,
    const AvalanchePassConfig& pc,
    std::vector<std::vector<double>>& signalOut,
    std::ofstream* pairsCsv, unsigned long printEvery,
    const std::string& tag,
    Component* scField, ComponentPoisson2d* scDeposit,
    const SpaceChargeConfig* scCfg,
    ChargeGrid* scGrid,
    std::vector<std::vector<double>>* promptSignalOut,
    std::vector<std::vector<double>>* delayedSignalOut) {
  const std::size_t nStrips = strips.size();
  bool useDelayedSignal = false;
  if (pc.enableSignal && pc.enableDelayedSignal) {
    // Probe a quarter of the way into the bulk so a silently-empty delayed
    // store is caught here rather than showing up as a missing slow lobe.
    const auto delayed = ValidateDelayedWeightingData(
        weightingCmp, strips, pc.yMin + 0.25 * (pc.yMax - pc.yMin));
    if (!delayed.available) {
      const std::string msg =
          "RunAvalanchePass: delayed signal requested, but " +
          delayed.message + ".";
      if (pc.requireDelayedWeightingData) {
        throw std::runtime_error(msg);
      }
      std::cerr << "WARNING: " << msg
                << " Continuing with prompt signal only.\n";
    } else {
      useDelayedSignal = true;
      if (!pc.delayedSignalTimesNs.empty()) {
        const double tMaxRequested =
            *std::max_element(pc.delayedSignalTimesNs.begin(),
                              pc.delayedSignalTimesNs.end());
        const double tEnd = pc.tStart + pc.tStep * pc.nBins;
        if (tMaxRequested > tEnd + 1.e-12) {
          std::cerr << "WARNING: delayed weighting data extend to "
                    << tMaxRequested << " ns, beyond the signal window ending at "
                    << tEnd << " ns; the slow lobe will be truncated.\n";
        }
      }
    }
  }
  signalOut.assign(nStrips, std::vector<double>(pc.nBins, 0.));
  if (promptSignalOut) {
    promptSignalOut->assign(nStrips, std::vector<double>(pc.nBins, 0.));
  }
  if (delayedSignalOut) {
    delayedSignalOut->assign(nStrips, std::vector<double>(pc.nBins, 0.));
  }
  unsigned long nTotal = 0;
  std::atomic<unsigned long> nDone{0};  // incremented lock-free by all threads
  const auto tStart = Clock::now();
  const long long nPrim = static_cast<long long>(primaries.size());

  // Garfield prints an oversubscription notice per avalanche when called
  // from inside an OpenMP region. That IS the intended behaviour here
  // (one Garfield thread per OpenMP thread), but it floods stderr.
  std::ofstream devNull;
  std::streambuf* oldCerr = nullptr;
  if (pc.silenceGarfieldStderr) {
    devNull.open("/dev/null");
    oldCerr = std::cerr.rdbuf(devNull.rdbuf());
  }

  // Garfield's TransportParticles progress meter is written to stdout.
  // Redirect stdout for the full pass so large MIP ensembles do not emit one
  // progress line per primary. main.cpp prints event summaries after return.
  std::ofstream devNullPass;
  std::streambuf* oldCoutPass = nullptr;
  if (pc.silenceGarfieldAvalancheStdout) {
    devNullPass.open("/dev/null");
    if (devNullPass) oldCoutPass = std::cout.rdbuf(devNullPass.rdbuf());
  }

  // Optional narrower setup-only suppression when full-pass suppression is
  // disabled.
  std::ofstream devNullSetup;
  std::streambuf* oldCout = nullptr;
  if (!oldCoutPass && pc.silenceGarfieldSetupStdout) {
    devNullSetup.open("/dev/null");
    if (devNullSetup) oldCout = std::cout.rdbuf(devNullSetup.rdbuf());
  }


  #pragma omp parallel if (!pc.forceSerial)
  {
    // Sensor construction prints to stdout; serialise it so 8 threads'
    // messages don't interleave into garbage. Once per thread, so free.
    Sensor localSensor;
    #pragma omp critical
    {
      localSensor.AddComponent(&driftCmp);
      // TCAD first: the sensor takes its medium from the first component
      // returning status 0. The screening perturbation is added second.
      if (scField) localSensor.AddComponent(scField);
      if (pc.enableSignal) {
        for (std::size_t k = 0; k < strips.size(); ++k) {
          // ComponentGrid ignores the label, so when a per-strip vector is
          // supplied each electrode must come from its own component.
          Component* wc = perStripWeighting.empty()
                              ? &weightingCmp
                              : perStripWeighting[k];
          localSensor.AddElectrode(wc, strips[k].label);
        }
        localSensor.SetTimeWindow(pc.tStart, pc.tStep, pc.nBins);
        if (useDelayedSignal) {
          localSensor.EnableDelayedSignal();
          if (!pc.delayedSignalTimesNs.empty()) {
            localSensor.SetDelayedSignalTimes(pc.delayedSignalTimesNs);
          }
          if (pc.delayedSignalAveragingOrder > 0) {
            localSensor.SetDelayedSignalAveragingOrder(
                pc.delayedSignalAveragingOrder);
          }
        }
      }
      localSensor.SetArea(pc.xMin, pc.yMin, pc.zMin,
                          pc.xMax, pc.yMax, pc.zMax);
    }

    // No worker may start avalanche processing until all local Sensors have
    // been built and stdout has been restored by one thread.
    #pragma omp barrier
    #pragma omp single
    {
      if (oldCout) {
        std::cout.flush();
        std::cout.rdbuf(oldCout);
        oldCout = nullptr;
      }
    }
    #pragma omp barrier

    AvalancheMC localAval;
    localAval.SetSensor(&localSensor);
    ConfigureAvalancheFast(localAval, pc.fineStepCm, pc.bulkStepCm,
                           pc.yFineLoCm, pc.yFineHiCm, pc.timeWindowNs,
                           pc.sizeCap, pc.enableSignal);
    localAval.EnableDiffusion(pc.diffusion);
    localAval.EnableMultiplication(pc.multiplication);
    // drift paths are stored only when we need them to deposit space
    // charge; they cost memory, so they stay off otherwise
    if ((scDeposit || scGrid) && scCfg) localAval.EnableDriftLines(true);

    // per-thread CSV buffer: writing inside the loop would take a lock on
    // every single pair, serialising all threads
    std::ostringstream localRows;

    /* Private grid per thread. ChargeGrid::Deposit writes to four nodes
       and is NOT thread-safe, so a shared grid would race; the copies are
       summed under one critical section at the end instead. */
    ChargeGrid localGrid;
    if (scGrid && scCfg) { localGrid = *scGrid; localGrid.Clear(); }

    #pragma omp for schedule(dynamic, pc.ompDynamicChunk) reduction(+ : nTotal)
    for (long long i = 0; i < nPrim; ++i) {
      const auto& p = primaries[i];
      localAval.AvalancheElectronHole(p[0], p[1], p[2], p[3]);
      // AddCharge is mutex-protected, so depositing from inside the
      // parallel loop is safe. Solve() is left to the caller.
      if (scDeposit && scCfg) {
        DepositAvalancheCharge(localAval, *scDeposit, *scCfg);
      }
      if (scGrid && scCfg) {
        DepositAvalancheIntoGrid(localAval, localGrid, *scCfg);
      }
      std::size_t ne = 0, ni = 0;
      localAval.GetAvalancheSize(ne, ni);
      nTotal += ne;
      if (pairsCsv) {
        localRows << p[0] * 1.e4 << "," << p[1] * 1.e4 << "," << ne << "\n";
      }
      if (printEvery) {
        const unsigned long done = ++nDone;  // atomic, no lock
        if (done % printEvery == 0) {
          #pragma omp critical
          {
            std::cout << "  [" << tag << "] pair " << done << "/" << nPrim
                      << "  [" << ElapsedS(tStart) << " s]" << std::endl;
          }
        }
      }
    }

    // superposition: summing each thread's induced signal is equivalent
    // to having accumulated them all in one sensor. Rows are flushed here
    // too, so mip_pairs.csv is grouped by thread rather than track order
    // (each row carries its own x,y so order carries no information).
    #pragma omp critical
    {
      if (pairsCsv) *pairsCsv << localRows.str();
      if (scGrid && scCfg) {
        for (std::size_t i = 0; i < scGrid->rho.size(); ++i) {
          scGrid->rho[i] += localGrid.rho[i];
        }
      }
      if (pc.enableSignal) {
        for (std::size_t k = 0; k < nStrips; ++k) {
          for (unsigned int b = 0; b < pc.nBins; ++b) {
            signalOut[k][b] += localSensor.GetSignal(strips[k].label, b);
            if (promptSignalOut) {
              (*promptSignalOut)[k][b] +=
                  localSensor.GetPromptSignal(strips[k].label, b);
            }
            if (delayedSignalOut) {
              (*delayedSignalOut)[k][b] +=
                  localSensor.GetDelayedSignal(strips[k].label, b);
            }
          }
        }
      }
    }
  }
  if (oldCout) std::cout.rdbuf(oldCout);
  if (oldCoutPass) {
    std::cout.flush();
    std::cout.rdbuf(oldCoutPass);
  }
  if (oldCerr) std::cerr.rdbuf(oldCerr);
  return nTotal;
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

/* impact-ionisation feedback comparison
- Runs explicit seed/feedback modes for the requested models.
- Uses adaptive statistics so heavy-tailed modes receive more events.
- Writes v2 CSV files rather than appending a new schema to legacy results. */
void RunModelComparison(
    AvalancheMC& avalLadder, MediumSilicon& si, double x0, double yInj,
    double stepCm, std::size_t ladderCap, const std::string& outDir,
    const std::string& biasLabel, double eMax, double gMax,
    const FeedbackScanConfig& cfg) {
  struct SeedResult {
    double meanNe = 0.;
    double semNe = 0.;
    double relativeSemNe = 0.;
    double medianNe = 0.;
    double meanNh = 0.;
    double semNh = 0.;
    double noiseFNe = 0.;
    std::size_t maxNe = 0;
    std::size_t maxNh = 0;
    int nCapped = 0;
    int nFailed = 0;
    int nDone = 0;
    bool capLimited = false;
    bool precisionReached = false;
    std::string status = "NOT_RUN";
    std::vector<std::size_t> neValues;
    std::vector<std::size_t> nhValues;
  };

  struct ModelResult {
    SeedResult eNoHoles;
    SeedResult eFull;
    SeedResult hFull;
    SeedResult ehFull;
  };

  if (cfg.models.empty()) {
    std::cerr << "RunModelComparison: no models requested.\n";
    return;
  }
  if (cfg.minEvents <= 0 || cfg.maxEvents < cfg.minEvents ||
      cfg.batchSize <= 0 || cfg.targetRelativeSem <= 0.) {
    std::cerr << "RunModelComparison: invalid FeedbackScanConfig.\n";
    return;
  }

  const auto setModel = [&](const std::string& model) -> bool {
    if (model == "massey") {
      si.SetImpactIonisationModelMassey();
    } else if (model == "okuto") {
      si.SetImpactIonisationModelOkutoCrowell();
    } else if (model == "grant") {
      si.SetImpactIonisationModelGrant();
    } else if (model == "vodm") {
      si.SetImpactIonisationModelVanOverstraetenDeMan();
    } else {
      std::cerr << "Unknown impact-ionisation model: " << model << "\n";
      return false;
    }
    return true;
  };

  const auto suffix = [&](const std::string& model,
                          const std::string& mode) {
    const std::string bias = biasLabel == "NA" ? "" : "_" + biasLabel + "V";
    return "feedback_" + model + "_" + mode + bias + ".txt";
  };

  const auto writeSizes = [&](const std::string& model,
                              const std::string& mode,
                              const std::vector<std::size_t>& values) {
    std::ofstream f(outDir + suffix(model, mode));
    for (const auto value : values) f << value << "\n";
  };

  const auto finalise = [&](SeedResult& r) {
    if (r.nDone <= 0) {
      r.status = "NO_EVENTS";
      return;
    }

    double sumNe = 0., sumNe2 = 0., sumNh = 0., sumNh2 = 0.;
    for (int i = 0; i < r.nDone; ++i) {
      const double ne = static_cast<double>(r.neValues[i]);
      const double nh = static_cast<double>(r.nhValues[i]);
      sumNe += ne;
      sumNe2 += ne * ne;
      sumNh += nh;
      sumNh2 += nh * nh;
    }

    r.meanNe = sumNe / r.nDone;
    r.meanNh = sumNh / r.nDone;
    const double varNe =
        std::max(0., sumNe2 / r.nDone - r.meanNe * r.meanNe);
    const double varNh =
        std::max(0., sumNh2 / r.nDone - r.meanNh * r.meanNh);
    r.semNe = std::sqrt(varNe / r.nDone);
    r.semNh = std::sqrt(varNh / r.nDone);
    r.relativeSemNe = r.meanNe > 0. ? r.semNe / r.meanNe : 0.;
    r.medianNe = MedianOf(r.neValues);
    r.noiseFNe = r.meanNe > 0.
        ? (sumNe2 / r.nDone) / (r.meanNe * r.meanNe) : 0.;
    r.precisionReached =
        r.meanNe == 0. || r.relativeSemNe <= cfg.targetRelativeSem;

    if (r.capLimited) {
      r.status = "CAP_LIMITED";
    } else if (!r.precisionReached) {
      r.status = r.noiseFNe >= cfg.heavyTailThresholdF
          ? "HEAVY_TAIL_UNRESOLVED" : "PRECISION_NOT_REACHED";
    } else if (r.noiseFNe >= cfg.heavyTailThresholdF) {
      r.status = "CONVERGED_HEAVY_TAIL";
    } else {
      r.status = "CONVERGED";
    }
  };

  const auto runSeedMode = [&](const std::string& label, auto&& launch) {
    SeedResult r;
    r.neValues.reserve(cfg.maxEvents);
    r.nhValues.reserve(cfg.maxEvents);
    const auto tStart = Clock::now();

    double sumNe = 0.;
    double sumNe2 = 0.;

    for (int i = 0; i < cfg.maxEvents; ++i) {
      const bool ok = launch(x0);
      if (!ok) ++r.nFailed;

      std::size_t ne = 0, nh = 0;
      avalLadder.GetAvalancheSize(ne, nh);
      r.neValues.push_back(ne);
      r.nhValues.push_back(nh);
      r.maxNe = std::max(r.maxNe, ne);
      r.maxNh = std::max(r.maxNh, nh);
      ++r.nDone;

      const double dne = static_cast<double>(ne);
      sumNe += dne;
      sumNe2 += dne * dne;

      if (ne >= ladderCap || nh >= ladderCap) ++r.nCapped;

      // Abort clearly cap-limited feedback branches before spending many
      // minutes sampling a quantity whose mean is set by the artificial cap.
      const bool hardTrip = r.nCapped >= 6 && r.nCapped * 5 > r.nDone;
      const bool softTrip = r.nDone >= 50 && r.nCapped * 10 > r.nDone;
      if (hardTrip || softTrip) {
        r.capLimited = true;
        std::cout << label << ": " << r.nCapped << "/" << r.nDone
                  << " events reached the avalanche cap; aborting mode.\n";
        break;
      }

      const bool checkpoint =
          r.nDone >= cfg.minEvents &&
          (r.nDone % cfg.batchSize == 0 || r.nDone == cfg.maxEvents);
      if (checkpoint) {
        const double mean = sumNe / r.nDone;
        const double var = std::max(0., sumNe2 / r.nDone - mean * mean);
        const double sem = std::sqrt(var / r.nDone);
        const double relSem = mean > 0. ? sem / mean : 0.;
        if (mean == 0. || relSem <= cfg.targetRelativeSem) break;
      }
    }

    finalise(r);
    std::cout << label
              << "   <Ne>=" << r.meanNe << " +- " << r.semNe
              << " (relSEM=" << 100. * r.relativeSemNe << "%, median="
              << r.medianNe << ")"
              << "   <Nh>=" << r.meanNh << " +- " << r.semNh
              << "   F_e=" << r.noiseFNe
              << "   max=(" << r.maxNe << "," << r.maxNh << ")"
              << "   capped=" << r.nCapped << "/" << r.nDone
              << "   failed=" << r.nFailed
              << "   status=" << r.status
              << "   [timer] " << ElapsedS(tStart) << " s\n";
    return r;
  };

  avalLadder.EnableMultiplication(true);
  std::vector<ModelResult> results(cfg.models.size());
  const auto tCmpStart = Clock::now();

  std::cout << "# feedback comparison: fine step = " << stepCm * 1.e7
            << " nm, cap = " << ladderCap
            << ", min/max events = " << cfg.minEvents << "/"
            << cfg.maxEvents
            << ", target relative SEM = " << 100. * cfg.targetRelativeSem
            << "%\n";
  std::cout << "# Note: in e_no_holes, generated holes are included in Nh "
               "but are not transported.\n";

  for (std::size_t im = 0; im < cfg.models.size(); ++im) {
    const std::string& model = cfg.models[im];
    if (!setModel(model)) continue;
    std::cout << "\n## model = " << model << "\n";

    ModelResult& r = results[im];
    r.eNoHoles = runSeedMode(model + " e_no_holes", [&](double xi) {
      return avalLadder.AvalancheElectron(xi, yInj, 0., 0., false);
    });
    r.eFull = runSeedMode(model + " e_full", [&](double xi) {
      return avalLadder.AvalancheElectron(xi, yInj, 0., 0., true);
    });
    if (cfg.runHoleSeed) {
      r.hFull = runSeedMode(model + " h_full", [&](double xi) {
        return avalLadder.AvalancheHole(xi, yInj, 0., 0., true);
      });
    }
    if (cfg.runPairSeed) {
      r.ehFull = runSeedMode(model + " eh_full", [&](double xi) {
        return avalLadder.AvalancheElectronHole(xi, yInj, 0., 0.);
      });
    }

    writeSizes(model, "e_no_holes_ne", r.eNoHoles.neValues);
    writeSizes(model, "e_full_ne", r.eFull.neValues);
    if (cfg.runHoleSeed) {
      writeSizes(model, "h_full_ne", r.hFull.neValues);
      writeSizes(model, "h_full_nh", r.hFull.nhValues);
    }
    if (cfg.runPairSeed) {
      writeSizes(model, "eh_full_ne", r.ehFull.neValues);
    }
  }

  const std::string detailPath = outDir + "feedback_results_v2.csv";
  std::ofstream detail(detailPath);
  detail << "bias,model,Epeak_line_Vcm,Epeak_global_Vcm,step_nm,cap,mode,"
            "N,meanNe,semNe,relativeSemNe,medianNe,meanNh,semNh,FNe,"
            "maxNe,maxNh,nCapped,nFailed,status\n";

  const auto writeDetail = [&](const std::string& model,
                               const std::string& mode,
                               const SeedResult& r) {
    detail << biasLabel << "," << model << "," << eMax << "," << gMax
           << "," << stepCm * 1.e7 << "," << ladderCap << "," << mode
           << "," << r.nDone << "," << r.meanNe << "," << r.semNe
           << "," << r.relativeSemNe << "," << r.medianNe << ","
           << r.meanNh << "," << r.semNh << "," << r.noiseFNe << ","
           << r.maxNe << "," << r.maxNh << "," << r.nCapped << ","
           << r.nFailed << "," << r.status << "\n";
  };

  const std::string summaryPath = outDir + "feedback_summary_v2.csv";
  std::ofstream summary(summaryPath);
  summary << "bias,model,eNoHoles,eFull,hFullNe,hFullNh,ehFull,"
             "secondaryFeedbackRatio,initialHoleRatio,eNoHolesStatus,"
             "eFullStatus,hFullStatus,ehFullStatus\n";

  std::cout << "\n# model  e_no_holes  e_full  h_full<Ne>  h_full<Nh>"
               "  eh_full  e_full/e_no_holes  eh_full/e_full  status\n";
  for (std::size_t im = 0; im < cfg.models.size(); ++im) {
    const auto& model = cfg.models[im];
    const auto& r = results[im];
    writeDetail(model, "e_no_holes", r.eNoHoles);
    writeDetail(model, "e_full", r.eFull);
    if (cfg.runHoleSeed) writeDetail(model, "h_full", r.hFull);
    if (cfg.runPairSeed) writeDetail(model, "eh_full", r.ehFull);

    const double secondaryFeedback = r.eNoHoles.meanNe > 0.
        ? r.eFull.meanNe / r.eNoHoles.meanNe : 0.;
    const double initialHoleRatio =
        cfg.runPairSeed && r.eFull.meanNe > 0.
        ? r.ehFull.meanNe / r.eFull.meanNe : 0.;
    const bool capLimited = r.eNoHoles.capLimited || r.eFull.capLimited ||
        (cfg.runHoleSeed && r.hFull.capLimited) ||
        (cfg.runPairSeed && r.ehFull.capLimited);

    std::cout << model << "  " << r.eNoHoles.meanNe << "  "
              << r.eFull.meanNe << "  " << r.hFull.meanNe << "  "
              << r.hFull.meanNh << "  " << r.ehFull.meanNe << "  "
              << secondaryFeedback << "  " << initialHoleRatio << "  "
              << (capLimited ? "CAP_LIMITED" : r.eFull.status) << "\n";

    summary << biasLabel << "," << model << "," << r.eNoHoles.meanNe
            << "," << r.eFull.meanNe << "," << r.hFull.meanNe << ","
            << r.hFull.meanNh << "," << r.ehFull.meanNe << ","
            << secondaryFeedback << "," << initialHoleRatio << ","
            << r.eNoHoles.status << "," << r.eFull.status << ","
            << r.hFull.status << "," << r.ehFull.status << "\n";
  }

  std::cout << "wrote " << detailPath << " and " << summaryPath << "\n";
  std::cout << "[timer] feedback comparison: " << ElapsedS(tCmpStart)
            << " s\n";
}

/* build the FEM region, grounded electrodes and mesh, then Initialise().
   Doping is deliberately left unset so this component carries only the
   deposited avalanche charge -- the TCAD map already holds the ionised
   dopants, and setting doping here would double-count them. */
bool SetupSpaceCharge(ComponentPoisson2d& sc, const SpaceChargeConfig& cfg,
                      Medium* medium) {
  if (cfg.xMaxCm <= cfg.xMinCm || cfg.yMaxCm <= cfg.yMinCm) {
    std::cerr << "SetupSpaceCharge: degenerate region." << std::endl;
    return false;
  }
  // rectangle: ComponentPoisson2d wants edges parallel to x or y,
  // anything else is approximated by a staircase
  const std::vector<double> xv = {cfg.xMinCm, cfg.xMaxCm,
                                  cfg.xMaxCm, cfg.xMinCm};
  const std::vector<double> yv = {cfg.yMinCm, cfg.yMinCm,
                                  cfg.yMaxCm, cfg.yMaxCm};
  if (!sc.AddRegion(xv, yv, medium, "bulk")) {
    std::cerr << "SetupSpaceCharge: AddRegion failed (does the medium have "
              << "a dielectric constant?)." << std::endl;
    return false;
  }
  // grounded Dirichlet boundaries. The electrode potentials themselves
  // live in the TCAD map; this component only carries the perturbation,
  // so both planes sit at 0 V.
  sc.AddElectrode(cfg.xMinCm, cfg.yMinCm, cfg.xMaxCm, cfg.yMinCm, 0., "top");
  sc.AddElectrode(cfg.xMinCm, cfg.yMaxCm, cfg.xMaxCm, cfg.yMaxCm, 0., "back");
  sc.SetRangeZ(cfg.zMinCm, cfg.zMaxCm);
  sc.SetMeshGranularity(cfg.hMinCm, cfg.hMaxCm);

  if (!sc.Initialise()) {
    std::cerr << "SetupSpaceCharge: Initialise() failed." << std::endl;
    return false;
  }
  // pre-solve the weighting potentials: the lazy path is mutex-guarded,
  // so leaving it to first use inside a parallel region would serialise
  // the threads on the first hit
  sc.RequestWeightingPotential("top");
  sc.RequestWeightingPotential("back");
  if (cfg.verbose) {
    std::cout << "  [spacecharge] Poisson mesh ready over x = ["
              << cfg.xMinCm * 1.e4 << ", " << cfg.xMaxCm * 1.e4 << "] um, y = ["
              << cfg.yMinCm * 1.e4 << ", " << cfg.yMaxCm * 1.e4 << "] um"
              << std::endl;
  }
  return true;
}

/* deposit one avalanche's carriers, weighted by residence time.
   AddCharge takes elementary charges per cm of depth, so the carrier
   count is divided by the effective z-extent. Electrons negative, holes
   positive: their spatial separation is what opposes the field. */
std::size_t DepositAvalancheCharge(const AvalancheMC& aval,
                                   ComponentPoisson2d& sc,
                                   const SpaceChargeConfig& cfg) {
  std::size_t nOutside = 0, nEmpty = 0;
  if (cfg.tWindowNs <= 0. || cfg.zExtentCm <= 0.) {
    std::cerr << "DepositAvalancheCharge: bad time window or z-extent."
              << std::endl;
    return 0;
  }
  const double norm = 1.0 / (cfg.zExtentCm * cfg.tWindowNs);

  auto deposit = [&](const std::vector<AvalancheMC::EndPoint>& eps,
                     const double sign) {
    for (const auto& ep : eps) {
      if (ep.path.size() < 2) { ++nEmpty; continue; }
      const double w = static_cast<double>(ep.weight) * sign * norm;
      for (std::size_t k = 0; k + 1 < ep.path.size(); ++k) {
        const auto& p0 = ep.path[k];
        const auto& p1 = ep.path[k + 1];
        const double dt = p1.t - p0.t;   // ns spent in this segment
        if (dt <= 0.) continue;
        if (!sc.AddCharge(0.5 * (p0.x + p1.x), 0.5 * (p0.y + p1.y), w * dt)) {
          ++nOutside;
        }
      }
    }
  };
  deposit(aval.GetElectrons(), -1.);
  deposit(aval.GetHoles(), +1.);

  if (nEmpty > 0 && cfg.verbose) {
    std::cerr << "  [spacecharge] WARNING: " << nEmpty << " carriers had no"
              << " stored path -- set aval.EnableDriftLines(true)"
              << std::endl;
  }
  return nOutside;
}

ScreeningFieldSample SampleScreeningField(
    ComponentPoisson2d& sc, const SpaceChargeConfig& cfg) {
  ScreeningFieldSample r;
  double ez = 0.;
  Medium* medium = nullptr;
  sc.ElectricField(cfg.xProbeCm, cfg.yProbeCm, 0., r.exVcm, r.eyVcm, ez,
                   r.potentialV, medium, r.status);
  r.magnitudeVcm = std::sqrt(r.exVcm * r.exVcm + r.eyVcm * r.eyVcm);
  r.depositedChargeEPerCm = sc.GetTotalDepositedCharge();
  return r;
}

ScreeningFieldGridSample SampleScreeningFieldGrid(
    ComponentPoisson2d& sc, const SpaceChargeConfig& cfg) {
  ScreeningFieldGridSample out;
  const int nx = std::max(1, cfg.fieldSampleNx);
  const int ny = std::max(1, cfg.fieldSampleNy);
  out.exVcm.reserve(static_cast<std::size_t>(nx) * ny);
  out.eyVcm.reserve(static_cast<std::size_t>(nx) * ny);

  const double xLo = std::max(cfg.xMinCm,
      cfg.xProbeCm - std::abs(cfg.fieldSampleXHalfWidthCm));
  const double xHi = std::min(cfg.xMaxCm,
      cfg.xProbeCm + std::abs(cfg.fieldSampleXHalfWidthCm));
  const double yLo = std::max(cfg.yMinCm,
      cfg.yProbeCm - std::abs(cfg.fieldSampleYHalfWidthCm));
  const double yHi = std::min(cfg.yMaxCm,
      cfg.yProbeCm + std::abs(cfg.fieldSampleYHalfWidthCm));

  double sum2 = 0.;
  for (int ix = 0; ix < nx; ++ix) {
    const double x = nx == 1 ? cfg.xProbeCm
        : xLo + (xHi - xLo) * ix / (nx - 1);
    for (int iy = 0; iy < ny; ++iy) {
      const double y = ny == 1 ? cfg.yProbeCm
          : yLo + (yHi - yLo) * iy / (ny - 1);
      double ex = 0., ey = 0., ez = 0., v = 0.;
      int status = 0;
      Medium* medium = nullptr;
      sc.ElectricField(x, y, 0., ex, ey, ez, v, medium, status);
      if (status != 0 || !std::isfinite(ex) || !std::isfinite(ey)) {
        out.exVcm.push_back(std::numeric_limits<double>::quiet_NaN());
        out.eyVcm.push_back(std::numeric_limits<double>::quiet_NaN());
        continue;
      }
      out.exVcm.push_back(ex);
      out.eyVcm.push_back(ey);
      ++out.nValid;
      const double mag2 = ex * ex + ey * ey;
      sum2 += mag2;
      out.maxMagnitudeVcm =
          std::max(out.maxMagnitudeVcm, std::sqrt(mag2));
    }
  }
  out.l2NormVcm = std::sqrt(sum2);
  return out;
}

ScreeningFieldChange CompareScreeningFieldGrids(
    const ScreeningFieldGridSample& previous,
    const ScreeningFieldGridSample& current) {
  ScreeningFieldChange out;
  if (previous.exVcm.size() != current.exVcm.size() ||
      previous.eyVcm.size() != current.eyVcm.size()) {
    return out;
  }

  double diff2 = 0., ref2 = 0.;
  double maxDiff = 0., maxRef = 0.;
  for (std::size_t i = 0; i < current.exVcm.size(); ++i) {
    const double pex = previous.exVcm[i];
    const double pey = previous.eyVcm[i];
    const double cex = current.exVcm[i];
    const double cey = current.eyVcm[i];
    if (!std::isfinite(pex) || !std::isfinite(pey) ||
        !std::isfinite(cex) || !std::isfinite(cey)) {
      continue;
    }
    const double dex = cex - pex;
    const double dey = cey - pey;
    const double d2 = dex * dex + dey * dey;
    const double r2 = pex * pex + pey * pey;
    diff2 += d2;
    ref2 += r2;
    maxDiff = std::max(maxDiff, std::sqrt(d2));
    maxRef = std::max(maxRef, std::sqrt(r2));
    ++out.nCompared;
  }
  if (out.nCompared == 0) return out;
  out.relativeL2 = std::sqrt(diff2) / std::max(1.e-12, std::sqrt(ref2));
  out.relativeMax = maxDiff / std::max(1.e-12, maxRef);
  return out;
}

/* screening-field diagnostic. A few thousand V/cm against a ~376 kV/cm
   gain layer is the expected scale. */
void ReportScreeningField(ComponentPoisson2d& sc,
                          const SpaceChargeConfig& cfg) {
  const auto r = SampleScreeningField(sc, cfg);
  std::cout << "  [spacecharge] deposited "
            << r.depositedChargeEPerCm << " e/cm; screening field at"
            << " probe = " << r.magnitudeVcm << " V/cm"
            << " (status " << r.status << ")" << std::endl;
}

/* ---- spatial charge mixing ------------------------------------------ */

void ChargeGrid::Init(const SpaceChargeConfig& cfg) {
  nx = std::max(2, cfg.mixNx);
  ny = std::max(2, cfg.mixNy);
  xMinCm = cfg.xMinCm; xMaxCm = cfg.xMaxCm;
  yMinCm = cfg.yMinCm; yMaxCm = cfg.yMaxCm;
  invDx = static_cast<double>(nx - 1) / (xMaxCm - xMinCm);
  invDy = static_cast<double>(ny - 1) / (yMaxCm - yMinCm);
  rho.assign(static_cast<std::size_t>(nx) * ny, 0.);
}

double ChargeGrid::NodeX(int ix) const {
  return xMinCm + (xMaxCm - xMinCm) * ix / (nx - 1);
}
double ChargeGrid::NodeY(int iy) const {
  return yMinCm + (yMaxCm - yMinCm) * iy / (ny - 1);
}
double ChargeGrid::Total() const {
  double t = 0.;
  for (double r : rho) t += r;
  return t;
}
double ChargeGrid::AbsoluteTotal() const {
  double t = 0.;
  for (double r : rho) t += std::abs(r);
  return t;
}

/* cloud-in-cell: split q among the four surrounding nodes by bilinear
   weight. Conserves charge exactly and keeps sub-cell position, unlike
   dumping the whole cell's charge at its centre. */
void ChargeGrid::Deposit(const double xCm, const double yCm,
                         const double q) {
  if (nx < 2 || ny < 2 || invDx <= 0. || invDy <= 0.) return;
  if (!Contains(xCm, yCm)) return;
  DepositUnchecked(xCm, yCm, q);
}

void ChargeGrid::DepositUnchecked(const double xCm, const double yCm,
                                  const double q) noexcept {
  const double gx = (xCm - xMinCm) * invDx;
  const double gy = (yCm - yMinCm) * invDy;
  const int ix = std::min(nx - 2, std::max(0, static_cast<int>(gx)));
  const int iy = std::min(ny - 2, std::max(0, static_cast<int>(gy)));
  const double fx = std::min(1., std::max(0., gx - ix));
  const double fy = std::min(1., std::max(0., gy - iy));

  const std::size_t i00 = static_cast<std::size_t>(iy) * nx + ix;
  const std::size_t i10 = i00 + 1;
  const std::size_t i01 = i00 + nx;
  const std::size_t i11 = i01 + 1;
  const double wx0 = 1. - fx;
  const double wy0 = 1. - fy;
  rho[i00] += q * wx0 * wy0;
  rho[i10] += q * fx  * wy0;
  rho[i01] += q * wx0 * fy;
  rho[i11] += q * fx  * fy;
}

std::size_t DepositAvalancheIntoGrid(const AvalancheMC& aval,
                                     ChargeGrid& grid,
                                     const SpaceChargeConfig& cfg) {
  std::size_t nOutside = 0, nEmpty = 0;
  if (cfg.tWindowNs <= 0. || cfg.zExtentCm <= 0.) return 0;
  const double norm = 1.0 / (cfg.zExtentCm * cfg.tWindowNs);

  auto deposit = [&](const std::vector<AvalancheMC::EndPoint>& eps,
                     const double sign) {
    for (const auto& ep : eps) {
      if (ep.path.size() < 2) { ++nEmpty; continue; }
      const double w = static_cast<double>(ep.weight) * sign * norm;
      for (std::size_t k = 0; k + 1 < ep.path.size(); ++k) {
        const auto& p0 = ep.path[k];
        const auto& p1 = ep.path[k + 1];
        const double dt = p1.t - p0.t;
        if (dt <= 0.) continue;
        const double xm = 0.5 * (p0.x + p1.x);
        const double ym = 0.5 * (p0.y + p1.y);
        if (!grid.Contains(xm, ym)) { ++nOutside; continue; }
        grid.DepositUnchecked(xm, ym, w * dt);
      }
    }
  };
  deposit(aval.GetElectrons(), -1.);
  deposit(aval.GetHoles(), +1.);
  if (nEmpty > 0 && cfg.verbose) {
    std::cerr << "  [spacecharge] WARNING: " << nEmpty << " carriers had no"
              << " stored path -- set aval.EnableDriftLines(true)\n";
  }
  return nOutside;
}

void LoadChargeGrid(ComponentPoisson2d& sc, const ChargeGrid& grid) {
  sc.ClearCharge();
  std::size_t nAccepted = 0, nRejected = 0;
  const double dx = grid.nx > 1
      ? (grid.xMaxCm - grid.xMinCm) / (grid.nx - 1) : 0.;
  const double dy = grid.ny > 1
      ? (grid.yMaxCm - grid.yMinCm) / (grid.ny - 1) : 0.;
  // Exact boundary coordinates may not belong to a FEM element under
  // strict lookup tests. Nudge boundary grid nodes a negligible distance
  // inside while retaining their charge and relative position.
  const double epsX = std::max(1.e-12, 1.e-6 * std::abs(dx));
  const double epsY = std::max(1.e-12, 1.e-6 * std::abs(dy));
  for (int iy = 0; iy < grid.ny; ++iy) {
    for (int ix = 0; ix < grid.nx; ++ix) {
      const double q = grid.rho[static_cast<std::size_t>(iy) * grid.nx + ix];
      if (q == 0.) continue;
      const double x = std::clamp(grid.NodeX(ix),
                                  grid.xMinCm + epsX,
                                  grid.xMaxCm - epsX);
      const double y = std::clamp(grid.NodeY(iy),
                                  grid.yMinCm + epsY,
                                  grid.yMaxCm - epsY);
      if (sc.AddCharge(x, y, q)) {
        ++nAccepted;
      } else {
        ++nRejected;
      }
    }
  }
  if (nAccepted == 0) {
    std::cerr << "  [spacecharge] WARNING: no charge accepted by the solver"
              << " -- is the mixing grid inside the FEM region?\n";
  } else if (nRejected > 0) {
    std::cerr << "  [spacecharge] WARNING: " << nRejected
              << " nonzero mixing-grid nodes were rejected by AddCharge.\n";
  }
}

void ApplyRelaxedCharge(ComponentPoisson2d& sc, ChargeGrid& mixed,
                        const ChargeGrid& fresh, const double lambda) {
  const double lam = std::min(1., std::max(0., lambda));
  if (mixed.rho.size() != fresh.rho.size()) {
    std::cerr << "ApplyRelaxedCharge: grid size mismatch.\n";
    return;
  }
  // Poisson is linear, so blending the density is exactly equivalent to
  // blending the resulting fields -- and unlike a scalar amplitude
  // rescale it is valid even when the cloud changes shape.
  for (std::size_t i = 0; i < mixed.rho.size(); ++i) {
    mixed.rho[i] = (1. - lam) * mixed.rho[i] + lam * fresh.rho[i];
  }
  LoadChargeGrid(sc, mixed);
}

