#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <random>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <memory>

#include "Garfield/ComponentUser.hh"
#include "Garfield/ComponentTcad2d.hh"
#include "Garfield/ComponentConstant.hh"
#include "Garfield/GeometrySimple.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/SolidBox.hh"
#include "Garfield/TrackHeed.hh"
#include "Garfield/TrackBichsel.hh"

#include "lgad_tools.hh"
#include "sheet_weighting.hh"

using namespace Garfield;


struct SignalMetrics {
  std::vector<double> integral;
  std::vector<double> positiveArea;
  std::vector<double> negativeArea;
  std::vector<double> absoluteArea;
  std::vector<double> peak;
  std::vector<double> positivePeak;
  std::vector<double> negativePeak;
  std::vector<double> peakTimeNs;
  std::vector<double> zeroCrossingAfterPeakNs;
  double totalIntegral = 0.;
  double totalPositiveArea = 0.;
  double totalNegativeArea = 0.;
  double totalAbsoluteArea = 0.;
  std::size_t peakStrip = 0;
};

static SignalMetrics MeasureSignal(
    const std::vector<std::vector<double>>& signal,
    const double tStartNs, const double tStepNs) {
  SignalMetrics r;
  const std::size_t nStrips = signal.size();
  r.integral.assign(nStrips, 0.);
  r.positiveArea.assign(nStrips, 0.);
  r.negativeArea.assign(nStrips, 0.);
  r.absoluteArea.assign(nStrips, 0.);
  r.peak.assign(nStrips, 0.);
  r.positivePeak.assign(nStrips, 0.);
  r.negativePeak.assign(nStrips, 0.);
  r.peakTimeNs.assign(nStrips, std::numeric_limits<double>::quiet_NaN());
  r.zeroCrossingAfterPeakNs.assign(
      nStrips, std::numeric_limits<double>::quiet_NaN());

  for (std::size_t k = 0; k < nStrips; ++k) {
    std::size_t iPeak = 0;
    for (std::size_t b = 0; b < signal[k].size(); ++b) {
      const double sample = signal[k][b];
      const double area = sample * tStepNs;
      r.integral[k] += area;
      r.absoluteArea[k] += std::abs(area);
      if (area > 0.) r.positiveArea[k] += area;
      if (area < 0.) r.negativeArea[k] += area;
      if (sample > r.positivePeak[k]) r.positivePeak[k] = sample;
      if (sample < r.negativePeak[k]) r.negativePeak[k] = sample;
      if (std::abs(sample) > std::abs(r.peak[k])) {
        r.peak[k] = sample;
        iPeak = b;
        r.peakTimeNs[k] = tStartNs + (b + 0.5) * tStepNs;
      }
    }

    if (!signal[k].empty() && r.peak[k] != 0.) {
      const double peakSign = std::copysign(1., r.peak[k]);
      for (std::size_t b = iPeak + 1; b < signal[k].size(); ++b) {
        const double y0 = signal[k][b - 1];
        const double y1 = signal[k][b];
        if (peakSign * y0 > 0. && peakSign * y1 <= 0.) {
          double fraction = 0.;
          const double dy = y1 - y0;
          if (dy != 0.) fraction = std::clamp(-y0 / dy, 0., 1.);
          r.zeroCrossingAfterPeakNs[k] =
              tStartNs + (static_cast<double>(b) - 0.5 + fraction) * tStepNs;
          break;
        }
      }
    }

    r.totalIntegral += r.integral[k];
    r.totalPositiveArea += r.positiveArea[k];
    r.totalNegativeArea += r.negativeArea[k];
    r.totalAbsoluteArea += r.absoluteArea[k];
    if (k > 0 && std::abs(r.peak[k]) > std::abs(r.peak[r.peakStrip])) {
      r.peakStrip = k;
    }
  }
  return r;
}


static std::string EventTag(const int eventId) {
  std::ostringstream out;
  out << "event" << std::setw(4) << std::setfill('0') << eventId;
  return out.str();
}

struct SampleStats {
  std::size_t n = 0;
  double mean = std::numeric_limits<double>::quiet_NaN();
  double sd = std::numeric_limits<double>::quiet_NaN();
  double sem = std::numeric_limits<double>::quiet_NaN();
};

static SampleStats ComputeSampleStats(const std::vector<double>& values) {
  SampleStats r;
  r.n = values.size();
  if (values.empty()) return r;
  double sum = 0.;
  for (const double x : values) sum += x;
  r.mean = sum / static_cast<double>(values.size());
  if (values.size() == 1) {
    r.sd = 0.;
    r.sem = 0.;
    return r;
  }
  double ss = 0.;
  for (const double x : values) {
    const double d = x - r.mean;
    ss += d * d;
  }
  r.sd = std::sqrt(ss / static_cast<double>(values.size() - 1));
  r.sem = r.sd / std::sqrt(static_cast<double>(values.size()));
  return r;
}


static void DumpCombinedFieldProfile(
    Component& tcad, ComponentPoisson2d& screening,
    const double xCm, const double yMinCm, const double yMaxCm,
    const int nPoints, const std::string& path,
    const bool announce = true) {
  std::ofstream out(path);
  if (!out) {
    std::cerr << "Could not open combined-field profile " << path << "\n";
    return;
  }
  out << "x_um,y_um,tcadEx_Vcm,tcadEy_Vcm,tcadEmag_Vcm,"
         "scEx_Vcm,scEy_Vcm,scEmag_Vcm,combinedEx_Vcm,combinedEy_Vcm,"
         "combinedEmag_Vcm,deltaEmag_Vcm,tcadStatus,scStatus\n";
  const int n = std::max(2, nPoints);
  for (int i = 0; i < n; ++i) {
    const double y = yMinCm + (yMaxCm - yMinCm) * i / (n - 1);
    double tx = 0., ty = 0., tz = 0., tv = 0.;
    double sx = 0., sy = 0., sz = 0., sv = 0.;
    Medium* tm = nullptr;
    Medium* sm = nullptr;
    int ts = 0, ss = 0;
    tcad.ElectricField(xCm, y, 0., tx, ty, tz, tv, tm, ts);
    screening.ElectricField(xCm, y, 0., sx, sy, sz, sv, sm, ss);
    const double tmag = std::sqrt(tx * tx + ty * ty);
    const double smag = std::sqrt(sx * sx + sy * sy);
    const double cx = tx + sx;
    const double cy = ty + sy;
    const double cmag = std::sqrt(cx * cx + cy * cy);
    out << xCm * 1.e4 << "," << y * 1.e4 << ","
        << tx << "," << ty << "," << tmag << ","
        << sx << "," << sy << "," << smag << ","
        << cx << "," << cy << "," << cmag << ","
        << (cmag - tmag) << "," << ts << "," << ss << "\n";
  }
  if (announce) std::cout << "wrote " << path << "\n";
}

struct RunConfig {
  // "diode" | "lgad" | "aclgad": sets .sta, readout layout and sheet model
  // together. Override with --device.
  std::string device = "aclgad";
  std::string file = "/home/ahaines561/HEP/MAS/Silvaco_dat/"
                        // "diode.sta";
                        // "lgad150V.sta";
                        // "lgad180V.sta";
                        // "lgad190V.sta";
                        // "aclgad230V.sta";
                          "aclgad300V.sta";
  std::string outDir = "output_files_aclgad_run/";

  double biasVOverride = std::numeric_limits<double>::quiet_NaN();
  double xTrackUm = 250.;
  bool doPositionScan = false;
  double scanXStartUm = 0.;
  double scanXEndUm = 500.;
  double scanXStepUm = 20.;
  int scanMipsPerPoint = 3;
  // Uniform random x instead of a grid (a grid aliases with the pad pitch).
  bool scanRandomUniform = false;
  int scanNRandomEvents = 200;
  unsigned scanRandomSeed = 20260819u;
  double scanCfdFraction = 0.2;
  bool scanUseGain = false;
  // Serial transport during the scan; set only if a crash recurs.
  bool scanForceSerial = false;
  bool doWeightingDump = false;
  // psi_k(x,y,t) maps, one file per time slice
  bool doDynamicWeightingDump = true;
  int weightingDumpNx = 200;
  int weightingDumpNy = 120;
  // Dump frame times; empty = production grid. Off-grid times are interpolated.
  std::vector<double> weightingDumpTimesNs = {};
  bool exitAfterWeightingDump = false;
  // Append DC contact + backplane columns (needed for sum_psi_residual).
  bool weightingDumpIncludeDcBackplane = true;
  bool doConvergenceScan = false;
  bool doFeedbackScan = false;
  bool doModelComparison = false;  // deprecated alias for doFeedbackScan
  bool doMIP = true;

  // Each gain event with delayed signal costs ~500 s.
  int nMips = 4;

  // Measure L_z from transport (replaces the fitted scZExtentUm).
  bool measureLz = true;
  // Gain-layer band for L_z; 0 = all depths.
  double lzBandHalfWidthUm = 2.5;
  int mipEventOffset = 0;
  bool mipRunStatic = true;
  bool mipRunScreened = false;
  bool mipWriteOverlaySignals = true;
  bool mipWritePerEventSignals = false;
  bool mipWritePairFiles = false;
  // Drift paths for tools/make_drift_gif.py (static pass, event 0 only).
  bool mipWriteDriftLines = false;
  long driftMaxCarriers = 300;
  bool mipWritePrimaries = false;

  int mipProgressEvery = 1;
  unsigned long mipPairProgressEvery = 0;
  int avalancheOmpChunk = 4;
  bool mipForceSerialAvalanches = false;
  bool EnableDiffusion = true;
  bool EnableMultiplication = true;
  // repeated avalanches on the frozen final field
  int mipFinalSampleCount = 0;
  bool mipFinalSampleEnableSignal = true;

  bool silenceGarfieldSensorSetup = true;

  // false: static analytic strips; true: imported TCAD dynamic weighting.
  bool useDynamicTcadWeighting = false;
  std::string weightingGridFile = "";
  std::string weightingBiasDataFile = "";
  double weightingDeltaV = 1.0;
  // true: voltage step dv*Theta(t); false: impulse-field maps dv*delta(t).
  bool dynamicWeightingUsesPotential = true;
  std::vector<double> delayedSignalTimesNs = {
      0., 0.02, 0.05, 0.1, 0.2, 0.5, 1., 2., 5., 10., 20., 50., 100.};
  // One vector per strip, aligned with delayedSignalTimesNs.
  std::vector<std::vector<std::string>> dynamicWeightingDataFiles;
  std::size_t delayedSignalAveragingOrder = 2;

  // Resistive-sheet (AC-coupled) weighting; gives the bipolar waveform.
  // false = prompt-only segmented-LGAD reference.
  bool useResistiveGridWeighting = true;
  // From aclgad.sta: N_s = 7.58e12 cm^-2 -> 2740 ohm/sq at mu_n = 300.
  // TODO: document the source of 1500.
  double sheetResistanceOhmSq = 1500.;
  double couplingEpsR = 3.9;  // SiO2, confirmed
  double couplingOxideThicknessUm = 0.1;
  // Ohmic DC contact: 0..3.0 um, measured from aclgad300V.sta.
  // Represented width is quantised to +-dx/2 (see tools/README.md).
  double dcContactX0Um = 0.;
  double dcContactX1Um = 3.0;
  // Grounded metal over oxide (rest of the `bulk` electrode): capacitive sink.
  std::vector<std::pair<double, double>> groundedPadsUm = {{3.0, 5.0}};
  std::size_t sheetNodes = 301;
  std::size_t sheetOutY = 81;
  std::size_t sheetStepsPerDecade = 60;
  /* Delayed weighting time grid.
     - Spacing must be >= signalStepNs: two slices in one bin over-count
       charge (Sensor.cc:859-868).
     - Placed greedily for <0.3% psi interpolation error; regenerate with
       scratchpad/psi_dump.cpp if R_sq, t_ox or pad layout change.
     - Ends at 200 ns: Delayed() clamps beyond the last slice.
     - Memory ~98 MB at 99 slices, nOutY=81, nSheet=301. */
  std::vector<double> resistiveWeightingTimesNs = {
      0.203, 0.416, 0.628, 0.856, 1.11, 1.36, 1.59, 1.85, 2.05, 2.27,
      2.52, 2.79, 3.1, 3.43, 3.81, 4.01, 4.22, 4.44, 4.67, 4.92,
      5.18, 5.45, 5.74, 6.05, 6.37, 6.7, 7.06, 7.43, 8.23, 9.13,
      9.61, 10.1, 10.6, 11.2, 11.8, 12.4, 13.1, 13.8, 14.5, 15.3,
      16.1, 16.9, 17.8, 18.8, 19.7, 20.8, 21.9, 23, 24.3, 25.5,
      26.9, 28.3, 29.8, 31.4, 33, 34.8, 36.6, 38.5, 40.6, 42.7,
      45, 47.4, 49.9, 52.5, 55.3, 58.2, 61.2, 67.9, 75.2, 87.8,
      108, 200
  };
  // Collect at the depletion edge; the n+ above it is undepleted.
  double depletionEdgeUm = 0.430466;

  double staticSignalWindowNs = 4.;
  // Must exceed the last weighting time plus drift time.
  double dynamicSignalWindowNs = 210.;
  // Delayed cost ~ 1/signalStepNs. Also the floor on weighting slice spacing.
  double signalStepNs = 0.2;

  std::string model = "okuto";
  double temperatureK = 300.0;
  /* "heed" or "bichsel". They agree to 0.3% on mean pairs
     (tools/ionisation_compare). Heed: matches arXiv:2505.05974, works in any
     material. Bichsel: Si only, per-electron positions, needs upstream
     TrackBichsel (2026-08-07). */
  std::string ionisationModel = "heed";

  std::vector<std::string> feedbackModels = {};  // empty = use cfg.model
  int feedbackMinEvents = 500;
  int feedbackMaxEvents = 5000;
  int feedbackBatchSize = 250;
  double feedbackTargetRelativeSem = 0.05;
  double feedbackHeavyTailThresholdF = 3.0;
  double feedbackTimeWindowNs = 6.;
  std::size_t feedbackSizeCap = 5000;
  bool feedbackRunHoleSeed = true;
  bool feedbackRunPairSeed = true;
  bool dumpTownsendTables = true;

  unsigned int siliconRegion = 0;
  bool assignAllRegionsToSilicon = false;  // diagnostic only
  double driftWindowNs = 4.;
  double fineStepNm = 75.;   // gain-layer band
  double bulkStepNm = 800.;  // elsewhere
  double fineBandHalfWidthUm = 2.5;
  double activeFieldMinVcm = 100.;

  bool spaceCharge = true;
  // Fallback L_z when useMeasuredLz is off. 0.21 was fitted; measured is ~5 um.
  double scZExtentUm = 0.21;
  // Requires measureLz and a static pass.
  bool useMeasuredLz = true;
  int scMaxIter = 16;
  double scRelaxation = 0.25;
  bool scWriteIterationHistory = false;
  // 0 quiet, 1 per iteration, 2 detailed
  int scConsoleVerbosity = 1;
  double scGainTol = 0.01;  // diagnostic only
  double scFieldTol = 0.01;
  int scStableIterations = 3;
  int scMixNx = 101;
  int scMixNy = 501;
  int scFieldSampleNx = 11;
  int scFieldSampleNy = 41;
  double scFieldSampleXHalfWidthUm = 5.;
  double scFieldSampleYHalfWidthUm = 1.;
  bool scWriteFinalFieldProfile = false;
  int scFinalFieldProfilePoints = 801;
  bool scFinalEvaluationPass = true;
  // Seed each event with the previous converged charge, scaled by pair count.
  bool scWarmStartEvents = true;
  double scWarmStartScaleMin = 0.5;
  double scWarmStartScaleMax = 2.0;
  bool scPrintStageTimings = true;
};

int main(int argc, char** argv) {
  RunConfig cfg;
  const auto tRunStart = Clock::now();

  // CLI covers run selection only; physics choices stay in RunConfig.
  // Dump frames are log-spaced over [tmin, tmax], clamped to the solved range.
  {
    double tMin = 0., tMax = 0.;
    long frames = 0;
    bool xTrackSetByUser = false;
    auto need = [&](int i, const char* what) {
      if (i + 1 >= argc) {
        std::cerr << what << " needs a value\n";
        std::exit(2);
      }
      return std::string(argv[i + 1]);
    };
    for (int i = 1; i < argc; ++i) {
      const std::string a = argv[i];
      if (a == "--help" || a == "-h") {
        std::cout <<
            "usage: main [--dump-only] [--dump-nx N] [--dump-ny N]\n"
            "            [--dump-frames N] [--dump-tmin ns] [--dump-tmax ns]\n"
            "            [--scan-random N] [--scan-gain] [--no-mip]\n"
            "            [--signal-step ns] [--device NAME] [--file PATH]\n"
            "  --device: diode | lgad | lgad_djte | aclgad (default)\n"
            "            [--x-track um] [--drift-lines] [--feedback-scan]\n"
            "            [--n-mips N]\n"
            "            [--out DIR]\n"
            "With no arguments, runs the configured job as before.\n";
        return 0;
      } else if (a == "--dump-only") {
        cfg.doDynamicWeightingDump = true;
        cfg.exitAfterWeightingDump = true;
      } else if (a == "--dump-nx") {
        cfg.weightingDumpNx = std::stoi(need(i, "--dump-nx")); ++i;
      } else if (a == "--dump-ny") {
        cfg.weightingDumpNy = std::stoi(need(i, "--dump-ny")); ++i;
      } else if (a == "--dump-frames") {
        frames = std::stol(need(i, "--dump-frames")); ++i;
      } else if (a == "--dump-tmin") {
        tMin = std::stod(need(i, "--dump-tmin")); ++i;
      } else if (a == "--dump-tmax") {
        tMax = std::stod(need(i, "--dump-tmax")); ++i;
      } else if (a == "--scan-random") {
        cfg.doPositionScan = true;
        cfg.scanRandomUniform = true;
        cfg.scanNRandomEvents = std::stoi(need(i, "--scan-random")); ++i;
      } else if (a == "--device") {
        cfg.device = need(i, "--device"); ++i;
        const std::string base = "/home/ahaines561/HEP/MAS/Silvaco_dat/";
        if (cfg.device == "diode") {
          cfg.file = base + "diode.sta";
          cfg.useResistiveGridWeighting = false;
        } else if (cfg.device == "lgad") {
          cfg.file = base + "lgad180V.sta";
          cfg.useResistiveGridWeighting = false;
        } else if (cfg.device == "lgad_djte") {
          cfg.file = base + "lgad_deepjte_180V_X80.sta";
          cfg.useResistiveGridWeighting = false;
        } else if (cfg.device == "aclgad") {
          cfg.file = base + "aclgad300V.sta";
          cfg.useResistiveGridWeighting = true;
        } else {
          std::cerr << "--device must be diode, lgad, lgad_djte or aclgad\n";
          return 2;
        }
      } else if (a == "--n-mips") {
        cfg.nMips = std::stoi(need(i, "--n-mips")); ++i;
      } else if (a == "--feedback-scan") {
        cfg.doFeedbackScan = true;
        cfg.doMIP = false;
      } else if (a == "--drift-lines") {
        cfg.mipWriteDriftLines = true;
        cfg.mipWritePairFiles = true;
      } else if (a == "--x-track") {
        cfg.xTrackUm = std::stod(need(i, "--x-track")); ++i;
        xTrackSetByUser = true;
      } else if (a == "--file") {
        cfg.file = need(i, "--file"); ++i;
      } else if (a == "--signal-step") {
        // must be <= the delayed grid's min spacing (0.2 ns)
        cfg.signalStepNs = std::stod(need(i, "--signal-step")); ++i;
      } else if (a == "--scan-gain") {
        cfg.scanUseGain = true;
      } else if (a == "--no-mip") {
        cfg.doMIP = false;
      } else if (a == "--out") {
        cfg.outDir = need(i, "--out"); ++i;
        if (cfg.outDir.back() != '/') cfg.outDir += '/';
      } else {
        std::cerr << "unknown argument '" << a << "' (try --help)\n";
        return 2;
      }
    }
    // On the 100 um devices the default lands in the anode/cathode gap; use
    // the anode centre.
    if (!xTrackSetByUser && cfg.device != "aclgad") cfg.xTrackUm = 23.5;
    if (frames > 1 || tMin > 0. || tMax > 0.) {
      const auto& prod = cfg.resistiveWeightingTimesNs;
      if (prod.empty()) {
        std::cerr << "--dump-frames/--dump-tmin/--dump-tmax need a solved "
                     "weighting grid; resistiveWeightingTimesNs is empty\n";
        return 2;
      }
      if (tMin <= 0.) tMin = prod.front();
      if (tMax <= 0.) tMax = prod.back();
      tMin = std::max(tMin, prod.front());
      tMax = std::min(tMax, prod.back());
      if (frames < 2) frames = static_cast<long>(prod.size());
      if (tMax <= tMin) {
        std::cerr << "--dump-tmax must exceed --dump-tmin\n";
        return 2;
      }
      cfg.weightingDumpTimesNs.clear();
      for (long i = 0; i < frames; ++i)
        cfg.weightingDumpTimesNs.push_back(
            tMin * std::pow(tMax / tMin, double(i) / double(frames - 1)));
      std::cout << "dump frames: " << frames << " log-spaced over " << tMin
                << " .. " << tMax << " ns\n";
    }
  }

  /* Two delayed slices in one signal bin over-count charge (Sensor.cc:859).
     Checked before the sheet solve so a mismatch fails fast. To regenerate:
       ./sheet_psi_dump 1500 0.1 > psi_dense.csv
       python3 ../tools/sheet_grid_design.py psi_dense.csv --tol 0.01 \
               --min-spacing <new signalStepNs> */
  if (cfg.useResistiveGridWeighting) {
    const auto& wt = cfg.resistiveWeightingTimesNs;
    if (wt.empty() || wt.front() < cfg.signalStepNs - 1.e-12) {
      std::cerr << "Delayed weighting grid must start at or after "
                   "signalStepNs (" << cfg.signalStepNs << " ns); first slice "
                   "is " << (wt.empty() ? 0. : wt.front()) << " ns.\n"
                   "  Either lower signalStepNs or regenerate the grid with "
                   "tools/sheet_grid_design.py --min-spacing "
                << cfg.signalStepNs << "\n";
      return 3;
    }
    double worst = 1e300;
    std::size_t iw = 0;
    for (std::size_t i = 1; i < wt.size(); ++i)
      if (wt[i] - wt[i - 1] < worst) { worst = wt[i] - wt[i - 1]; iw = i; }
    if (worst < cfg.signalStepNs - 1.e-12) {
      std::cerr << "Delayed weighting times " << wt[iw - 1] << " and " << wt[iw]
                << " ns are " << worst << " ns apart, closer than one signal "
                   "bin (" << cfg.signalStepNs << " ns). Sensor would "
                   "over-count the delayed charge.\n"
                   "  Either set signalStepNs <= " << worst
                << " or regenerate the grid with tools/sheet_grid_design.py "
                   "--min-spacing " << cfg.signalStepNs << "\n";
      return 3;
    }
  }

  double biasV = cfg.biasVOverride;
  if (std::isnan(biasV)) biasV = ParseBiasFromFilename(cfg.file);
  const std::string biasLabel = FormatBias(biasV);

  std::filesystem::create_directories(cfg.outDir);
  std::cout << "bias = " << biasLabel << (std::isnan(biasV) ? "" : " V")
            << ", output directory = " << cfg.outDir << std::endl;

  ComponentTcad2d cmp;
  if (!cmp.InitialiseSilvaco(cfg.file)) return 1;

  MediumSilicon si;
  si.SetTemperature(cfg.temperatureK);
  if (!SetImpactIonisationModel(si, cfg.model)) return 2;
  std::cout << "impact-ionisation model: " << cfg.model
            << ", temperature = " << cfg.temperatureK << " K"
            << std::endl;

  if (cmp.GetNumberOfRegions() == 0) {
    std::cerr << "Field map contains no regions.\n";
    return 1;
  }
  if (cfg.assignAllRegionsToSilicon) {
    std::cout << "WARNING: assigning every imported region to silicon; "
                 "use only for a controlled diagnostic.\n";
    for (unsigned int i = 0; i < cmp.GetNumberOfRegions(); ++i) {
      cmp.SetMedium(i, &si);
    }
  } else {
    if (cfg.siliconRegion >= cmp.GetNumberOfRegions()) {
      std::cerr << "Requested silicon region " << cfg.siliconRegion
                << " but the map has only " << cmp.GetNumberOfRegions()
                << " regions.\n";
      return 1;
    }
    cmp.SetMedium(cfg.siliconRegion, &si);
  }
  cmp.SetRangeZ(-5.e-4, 5.e-4);
  cmp.PrintRegions();

  double bx0 = 0., by0 = 0., bz0 = 0., bx1 = 0., by1 = 0., bz1 = 0.;
  if (!cmp.GetBoundingBox(bx0, by0, bz0, bx1, by1, bz1)) {
    std::cerr << "Could not get the bounding box from the field map;\n"
              << "set the device extent manually in the macro.\n";
    return 1;
  }
  std::cout << "map extent: x = [" << bx0 * 1.e4 << ", " << bx1 * 1.e4
            << "] um, y = [" << by0 * 1.e4 << ", " << by1 * 1.e4
            << "] um" << std::endl;

  ScanValidity(cmp, bx0, bx1, by0, by1, 400, 400,
              cfg.outDir + "validity_map.csv");

  double x0 = cfg.xTrackUm * 1.e-4 + 0.13e-4;
  if (x0 <= bx0 || x0 >= bx1) {
    std::cout << "requested x = " << cfg.xTrackUm
              << " um is outside the map; using the device centre."
              << std::endl;
    x0 = 0.5 * (bx0 + bx1) + 0.13e-4;
  }
  std::cout << "probe/track line at x = " << x0 * 1.e4 << " um"
            << std::endl;

  const auto prof = ScanFieldProfile(cmp, x0, by0, by1, 400,
                                     cfg.outDir + (biasLabel == "NA"
                                         ? "profile.csv"
                                         : "profile_" + biasLabel + "V.csv"),
                                     cfg.activeFieldMinVcm);
  if (!prof.valid) return 1;
  const double yTop = prof.yTop;
  // depletion edge, from the .sta electron density (aclgad_extract.py)
  const double yCollect = std::max(yTop + 0.02e-4,
                                   cfg.depletionEdgeUm * 1.e-4);
  const double yBot = prof.yBot;
  const double d = prof.d;
  const double yGain = prof.yGain;
  const double eMax = prof.eMax;

  // 2D hotspot search, seeded with the probe line's peak
  double globalMaxField = eMax;
  double globalMaxX = x0;
  double globalMaxY = yGain;
  const auto tHotspotStart = Clock::now();
  for (int ix = 0; ix <= 60; ++ix) {
    const double x = bx0 + ((bx1 - bx0) * ix) / 60. + 0.077e-4;
    for (int iy = 0; iy <= 150; ++iy) {
      const double y = by0 + ((by1 - by0) * iy) / 150. + 0.0113e-4;
      double ex, ey, ez, v; int st; Medium* m = nullptr;
      cmp.ElectricField(x, y, 0., ex, ey, ez, v, m, st);
      if (st != 0) continue;
      const double e = std::sqrt(ex * ex + ey * ey);
      if (e > globalMaxField) { globalMaxField = e; globalMaxX = x; globalMaxY = y; }
    }
    for (int iy = 0; iy <= 120; ++iy) {
      const double y = yTop - 0.5e-4 + (3.5e-4 * iy) / 120.;
      double ex, ey, ez, v; int st; Medium* m = nullptr;
      cmp.ElectricField(x, y, 0., ex, ey, ez, v, m, st);
      if (st != 0) continue;
      const double e = std::sqrt(ex * ex + ey * ey);
      if (e > globalMaxField) { globalMaxField = e; globalMaxX = x; globalMaxY = y; }
    }
  }
  std::cout << "global max field " << globalMaxField << " V/cm at (x, y) = ("
            << globalMaxX * 1.e4 << ", " << globalMaxY * 1.e4 << ") um"
            << std::endl;
  std::cout << "[timer] hotspot scan: " << ElapsedS(tHotspotStart)
            << " s" << std::endl;

  // x = xTrack +- 5 um, full depth
  {
    const auto tDumpStart = Clock::now();
    const double xLo = std::max(bx0, (cfg.xTrackUm - 5.) * 1.e-4);
    const double xHi = std::min(bx1, (cfg.xTrackUm + 5.) * 1.e-4);
    const double dxCm = 0.1e-4, dyCm = 0.05e-4;
    const int nx = static_cast<int>(std::round((xHi - xLo) / dxCm));
    const int ny = static_cast<int>(std::round((by1 - by0) / dyCm));
    const auto r = DumpElectricField(cmp, xLo, xHi, by0, by1, nx, ny,
                                     cfg.outDir + "efield_window.csv");
    std::cout << "field dump: x = [" << xLo * 1.e4 << ", " << xHi * 1.e4
              << "] um -> efield_window.csv (" << r.nRows << " points, "
              << r.nBad << " invalid skipped)" << std::endl;
    std::cout << "  window max |E| = " << r.wMax << " V/cm at (x, y) = ("
              << r.wMaxX * 1.e4 << ", " << r.wMaxY * 1.e4 << ") um"
              << std::endl;
    std::cout << "[timer] field dump: " << ElapsedS(tDumpStart) << " s"
              << std::endl;
  }
  {
    const auto tFullStart = Clock::now();
    const auto r = DumpElectricField(cmp, bx0, bx1, by0, by1, 250, 250,
                                     cfg.outDir + "efield_full.csv");
    std::cout << "full-device dump -> efield_full.csv (" << r.nRows
              << " points, " << r.nBad << " invalid skipped)" << std::endl;
    std::cout << "[timer] full-device dump: " << ElapsedS(tFullStart)
              << " s" << std::endl;
  }

  // diode/LGAD: anode 0..47, cathode 53..100 um (from .sta)
  const std::vector<ReadoutStrip> strips =
      (cfg.device == "aclgad")
          ? std::vector<ReadoutStrip>{{"strip0", 50., 20.},
                                      {"strip1", 245., 25.},
                                      {"strip2", 450., 20.}}
          : std::vector<ReadoutStrip>{{"anode", 23.5, 23.5},
                                      {"cathode", 76.5, 23.5}};
  std::cout << "device = " << cfg.device << ", " << strips.size()
            << " readout electrodes, resistive weighting "
            << (cfg.useResistiveGridWeighting ? "ON" : "OFF") << "\n";
  if (!StripsInsideMap(strips, bx0, bx1)) return 1;

  ComponentAnalyticField analyticWeighting;
  ComponentTcad2d dynamicWeighting;
  ComponentUser resistiveWeighting;
  std::unique_ptr<mass::SheetWeighting> sheetModel;
  Component* weightingCmp = nullptr;
  std::vector<Component*> perStripWeighting;
  bool dynamicWeightingActive = false;

  if (cfg.useDynamicTcadWeighting) {
    if (cfg.weightingGridFile.empty() ||
        cfg.weightingBiasDataFile.empty()) {
      std::cerr << "Dynamic weighting requested, but weightingGridFile or "
                   "weightingBiasDataFile is empty.\n";
      return 3;
    }
    if (cfg.dynamicWeightingDataFiles.size() != strips.size()) {
      std::cerr << "Dynamic weighting requires one time-slice file vector "
                   "per readout strip.\n";
      return 3;
    }
    if (!dynamicWeighting.Initialise(cfg.weightingGridFile,
                                     cfg.weightingBiasDataFile)) {
      std::cerr << "Could not initialise the TCAD weighting component.\n";
      return 3;
    }

    for (std::size_t k = 0; k < strips.size(); ++k) {
      if (cfg.dynamicWeightingDataFiles[k].size() !=
          cfg.delayedSignalTimesNs.size()) {
        std::cerr << "Dynamic weighting file count for " << strips[k].label
                  << " does not match delayedSignalTimesNs.\n";
        return 3;
      }
      for (std::size_t it = 0; it < cfg.delayedSignalTimesNs.size(); ++it) {
        const double tNs = cfg.delayedSignalTimesNs[it];
        const std::string& dataFile =
            cfg.dynamicWeightingDataFiles[k][it];
        const bool ok = cfg.dynamicWeightingUsesPotential
            ? dynamicWeighting.SetDynamicWeightingPotential(
                  cfg.weightingBiasDataFile, dataFile,
                  cfg.weightingDeltaV, tNs, strips[k].label)
            : dynamicWeighting.SetDynamicWeightingField(
                  cfg.weightingBiasDataFile, dataFile,
                  cfg.weightingDeltaV, tNs, strips[k].label);
        if (!ok) {
          std::cerr << "Failed to load dynamic weighting slice for "
                    << strips[k].label << " at t=" << tNs
                    << " ns from " << dataFile << "\n";
          return 3;
        }
      }
    }

    const double wProbeY = yTop + 0.25 * (yBot - yTop);
    const auto validation =
        ValidateDelayedWeightingData(dynamicWeighting, strips, wProbeY);
    if (!validation.available) {
      std::cerr << "Dynamic weighting validation failed: "
                << validation.message << "\n";
      return 3;
    }
    weightingCmp = &dynamicWeighting;
    dynamicWeightingActive = true;
    std::cout << "dynamic TCAD weighting enabled with "
              << validation.timesNs.size() << " time slices through "
              << validation.timesNs.back() << " ns\n";
  } else if (cfg.useResistiveGridWeighting) {
    mass::SheetConfig sc;
    sc.xMin = bx0;
    sc.xMax = bx1;
    sc.ySheet = cfg.depletionEdgeUm * 1.e-4;
    sc.yBack = yBot;
    sc.tOx = cfg.couplingOxideThicknessUm * 1.e-4;
    sc.epsOx = cfg.couplingEpsR;
    sc.rSheet = cfg.sheetResistanceOhmSq;
    sc.dcX0 = cfg.dcContactX0Um * 1.e-4;
    sc.dcX1 = cfg.dcContactX1Um * 1.e-4;
    sc.groundedPads.clear();
    for (const auto& g : cfg.groundedPadsUm)
      sc.groundedPads.emplace_back(g.first * 1.e-4, g.second * 1.e-4);
    sc.nSheet = cfg.sheetNodes;
    sc.nOutY = cfg.sheetOutY;
    sc.stepsPerDecade = cfg.sheetStepsPerDecade;
    sc.timesNs = cfg.resistiveWeightingTimesNs;
    sc.pads.clear();
    for (const auto& st : strips) {
      sc.pads.emplace_back((st.centerUm - st.halfWidthUm) * 1.e-4,
                           (st.centerUm + st.halfWidthUm) * 1.e-4);
    }
    std::cout << "solving resistive-sheet weighting: R_sq = " << sc.rSheet
              << " ohm/sq, eps_ox = " << sc.epsOx << ", sheet at y = "
              << cfg.depletionEdgeUm << " um, " << sc.nSheet << " nodes, "
              << sc.timesNs.size() << " time slices\n";
    sheetModel = std::make_unique<mass::SheetWeighting>(sc);
    if (!sheetModel->SelfTest(std::cout)) {
      std::cerr << "resistive-sheet self test FAILED\n";
      return 3;
    }

    // model is immutable, so these lambdas are thread-safe
    const mass::SheetWeighting* model = sheetModel.get();
    for (std::size_t k = 0; k < strips.size(); ++k) {
      resistiveWeighting.SetWeightingPotential(
          [model, k](const double x, const double y, const double) {
            return model->Prompt(k, x, y);
          },
          strips[k].label);
      resistiveWeighting.SetDelayedWeightingPotential(
          [model, k](const double x, const double y, const double,
                     const double t) {
            return model->Delayed(k, x, y, t);
          },
          strips[k].label);
    }
    // DC contact + backplane: dump only, not readout channels
    if (cfg.weightingDumpIncludeDcBackplane) {
      const std::size_t kDC = strips.size(), kBP = strips.size() + 1;
      for (auto kv : {std::make_pair(kDC, std::string("dc")),
                      std::make_pair(kBP, std::string("backplane"))}) {
        const std::size_t k = kv.first;
        resistiveWeighting.SetWeightingPotential(
            [model, k](const double x, const double y, const double) {
              return model->Prompt(k, x, y);
            },
            kv.second);
        resistiveWeighting.SetDelayedWeightingPotential(
            [model, k](const double x, const double y, const double,
                       const double t) {
              return model->Delayed(k, x, y, t);
            },
            kv.second);
      }
    }
    {
    }
    resistiveWeighting.SetDelayedSignalTimes(cfg.resistiveWeightingTimesNs);
    weightingCmp = &resistiveWeighting;
    dynamicWeightingActive = true;

    const double wProbeY = yTop + 0.25 * (yBot - yTop);
    const auto v =
        ValidateDelayedWeightingData(resistiveWeighting, strips, wProbeY);
    if (!v.available) {
      std::cerr << "resistive weighting validation failed: " << v.message
                << "\n";
      return 3;
    }
    std::cout << "resistive-sheet AC-coupled weighting enabled: "
              << v.timesNs.size() << " delayed slices through "
              << v.timesNs.back() << " ns\n";
  } else {
    analyticWeighting.AddPlaneY(yTop, 1., "top");
    analyticWeighting.AddPlaneY(yBot, 0., "back");
    for (const auto& s : strips) {
      const double xc = s.centerUm * 1.e-4;
      const double hw = s.halfWidthUm * 1.e-4;
      analyticWeighting.AddStripOnPlaneY(
          'z', yTop, xc - hw, xc + hw, s.label);
    }
    weightingCmp = &analyticWeighting;
    std::cout << "static analytic weighting enabled: prompt-only LGAD "
                 "reference, not a dynamic AC-LGAD waveform\n";
  }

  PrintWeightingSanity(*weightingCmp, strips, yTop, yCollect, yBot);
  if (cfg.doDynamicWeightingDump && dynamicWeightingActive) {
    const auto& dumpTimes = cfg.weightingDumpTimesNs.empty()
                                ? cfg.resistiveWeightingTimesNs
                                : cfg.weightingDumpTimesNs;
    DumpDynamicWeighting(*weightingCmp, strips, bx0, bx1, yCollect, yBot,
                         dumpTimes,
                         cfg.outDir + "weighting" +
                             (biasLabel == "NA" ? "" : "_" + biasLabel + "V"),
                         cfg.weightingDumpNx, cfg.weightingDumpNy,
                         cfg.weightingDumpIncludeDcBackplane
                             ? std::vector<std::string>{"dc", "backplane"}
                             : std::vector<std::string>{});
    if (cfg.exitAfterWeightingDump) {
      std::cout << "exitAfterWeightingDump set: stopping before transport\n";
      return 0;
    }
  }

  const double signalWindowNs = dynamicWeightingActive
      ? cfg.dynamicSignalWindowNs : cfg.staticSignalWindowNs;
  if (cfg.signalStepNs <= 0. || signalWindowNs <= 0.) {
    std::cerr << "Invalid signal time step or signal window.\n";
    return 3;
  }
  const unsigned int signalBins = static_cast<unsigned int>(
      std::ceil(signalWindowNs / cfg.signalStepNs));
  const std::vector<double>& activeWeightingTimes =
      cfg.useResistiveGridWeighting ? cfg.resistiveWeightingTimesNs
                                    : cfg.delayedSignalTimesNs;
  if (dynamicWeightingActive &&
      !activeWeightingTimes.empty() &&
      signalWindowNs + 1.e-12 < activeWeightingTimes.back()) {
    std::cerr << "Dynamic signal window (" << signalWindowNs
              << " ns) is shorter than the last weighting-map time ("
              << activeWeightingTimes.back() << " ns).\n";
    return 3;
  }

  Sensor sensor;
  {
    std::ofstream nullOut;
    std::streambuf* oldCout = nullptr;
    if (cfg.silenceGarfieldSensorSetup) {
      nullOut.open("/dev/null");
      if (nullOut) oldCout = std::cout.rdbuf(nullOut.rdbuf());
    }
    sensor.AddComponent(&cmp);
    for (std::size_t k = 0; k < strips.size(); ++k) {
      Component* wc = perStripWeighting.empty() ? weightingCmp
                                                : perStripWeighting[k];
      sensor.AddElectrode(wc, strips[k].label);
    }
    sensor.SetTimeWindow(0., cfg.signalStepNs, signalBins);
    if (dynamicWeightingActive) {
      sensor.EnableDelayedSignal();
      sensor.SetDelayedSignalTimes(activeWeightingTimes);
      sensor.SetDelayedSignalAveragingOrder(
          cfg.delayedSignalAveragingOrder);
    }
    sensor.SetArea(bx0, yCollect, -5.e-4, bx1, yBot, 5.e-4);
    if (oldCout) {
      std::cout.flush();
      std::cout.rdbuf(oldCout);
    }
  }
  std::cout << "Sensor has " << sensor.GetNumberOfElectrodes()
            << " electrodes (expect " << strips.size() << ")" << std::endl;

  if (cfg.doWeightingDump) {
    DumpWeightingField(*weightingCmp, strips, bx0, bx1, yTop, yBot,
                       cfg.outDir + "wfield_full.txt");
  }

  const double fineStepCm = cfg.fineStepNm * 1.e-7;
  const double bulkStepCm = cfg.bulkStepNm * 1.e-7;
  const double fineBandCm = cfg.fineBandHalfWidthUm * 1.e-4;
  const double yFineLo = std::max(yGain - fineBandCm, yTop);
  const double yFineHi = std::min(yGain + fineBandCm, yBot);

  const std::size_t sizeCap = 5000;

  const double yInj = std::min(yGain + 5.e-4, 0.5 * (yGain + yBot));
  double ladderStepCm = fineStepCm;
  std::atomic<long long> ladderCalls{0}, ladderFine{0}, ladderCoarse{0},
      ladderPrinted{0};
  AvalancheMC avalLadder;
  avalLadder.SetSensor(&sensor);
  const bool controlledFeedbackScan =
      cfg.doFeedbackScan || cfg.doModelComparison;
  const std::size_t ladderCap = controlledFeedbackScan
      ? cfg.feedbackSizeCap : sizeCap;
  const double ladderTimeWindowNs = controlledFeedbackScan
      ? cfg.feedbackTimeWindowNs : cfg.driftWindowNs;
  ConfigureAvalanche(avalLadder, ladderStepCm, bulkStepCm, yFineLo, yFineHi,
                     ladderTimeWindowNs, ladderCap, /*enableSignal=*/false,
                     /*multithreading=*/!controlledFeedbackScan,
                     controlledFeedbackScan ? "feedback" : "ladder",
                     ladderCalls, ladderFine, ladderCoarse, ladderPrinted);

  if (cfg.doConvergenceScan) {
    RunConvergenceScan(avalLadder, ladderStepCm, x0, yInj, ladderCap,
                       cfg.outDir);
    std::cout << "[stepfn ladder] calls=" << ladderCalls.load()
              << " fine=" << ladderFine.load() << " coarse="
              << ladderCoarse.load() << std::endl;
  }

  if (cfg.doFeedbackScan || cfg.doModelComparison) {
    if (cfg.doModelComparison && !cfg.doFeedbackScan) {
      std::cout << "NOTE: doModelComparison is a deprecated alias; "
                   "running the explicit feedback scan.\n";
    }

    FeedbackScanConfig feedbackCfg;
    feedbackCfg.models = cfg.feedbackModels.empty()
        ? std::vector<std::string>{cfg.model} : cfg.feedbackModels;
    feedbackCfg.minEvents = cfg.feedbackMinEvents;
    feedbackCfg.maxEvents = cfg.feedbackMaxEvents;
    feedbackCfg.batchSize = cfg.feedbackBatchSize;
    feedbackCfg.targetRelativeSem = cfg.feedbackTargetRelativeSem;
    feedbackCfg.heavyTailThresholdF =
        cfg.feedbackHeavyTailThresholdF;
    feedbackCfg.runHoleSeed = cfg.feedbackRunHoleSeed;
    feedbackCfg.runPairSeed = cfg.feedbackRunPairSeed;

    if (cfg.dumpTownsendTables) {
      for (const auto& feedbackModel : feedbackCfg.models) {
        if (!SetImpactIonisationModel(si, feedbackModel)) return 2;
        DumpTownsendCoefficients(
            si, cfg.outDir + "townsend_garfield_" + feedbackModel +
                    "_" + std::to_string(static_cast<int>(
                        std::lround(cfg.temperatureK))) + "K.csv");
      }
    }

    RunModelComparison(avalLadder, si, x0, yInj, ladderStepCm, ladderCap,
                       cfg.outDir, biasLabel, eMax, globalMaxField,
                       feedbackCfg);
    std::cout << "[stepfn feedback] calls=" << ladderCalls.load()
              << " fine=" << ladderFine.load() << " coarse="
              << ladderCoarse.load() << std::endl;

    if (!SetImpactIonisationModel(si, cfg.model)) return 2;
  }

  // The position scan reuses the MIP setup below, so don't return early for it.
  if (!cfg.doMIP && !cfg.doPositionScan) {
    std::cout << "doMIP=false -- diagnostics complete; skipping the MIP pass."
              << std::endl;
    return 0;
  }
  const int nMipRun = cfg.doMIP ? cfg.nMips : 0;
  if (!cfg.doMIP)
    std::cout << "doMIP=false -- running the position scan only.\n";

  if (cfg.doMIP && cfg.nMips <= 0) {
    std::cerr << "nMips must be positive.\n";
    return 2;
  }
  const bool runStatic = cfg.mipRunStatic;
  const bool runScreened = cfg.mipRunScreened && cfg.spaceCharge;
  if (cfg.mipRunScreened && !cfg.spaceCharge) {
    std::cout << "NOTE: mipRunScreened=true but spaceCharge=false; "
                 "screened mode is disabled.\n";
  }
  if (!runStatic && !runScreened) {
    std::cerr << "No MIP mode selected. Enable mipRunStatic and/or "
                 "mipRunScreened.\n";
    return 2;
  }

  SolidBox box(0.5 * (bx0 + bx1), 0.5 * (yTop + yBot), 0.,
               0.5 * (bx1 - bx0), 0.5 * d, 5.e-4);
  GeometrySimple geo;
  geo.AddSolid(&box, &si);
  ComponentConstant cmpHeed;
  cmpHeed.SetGeometry(&geo);
  cmpHeed.SetElectricField(0., 100., 0.);
  Sensor heedSensor;
  heedSensor.AddComponent(&cmpHeed);

  TrackHeed track;
  track.SetSensor(&heedSensor);
  track.SetParticle("pi");
  track.SetMomentum(180.e9);

  const bool useBichsel = (cfg.ionisationModel == "bichsel");
  if (cfg.ionisationModel != "heed" && !useBichsel) {
    std::cerr << "ionisationModel must be \"heed\" or \"bichsel\", got \""
              << cfg.ionisationModel << "\"\n";
    return 1;
  }
  // Lazily built: the constructor throws if GARFIELD_INSTALL is unset.
  std::unique_ptr<TrackBichsel> bichselOwner;
  if (useBichsel) {
    try {
      bichselOwner = std::make_unique<TrackBichsel>(&heedSensor);
    } catch (const std::exception& e) {
      std::cerr << "TrackBichsel could not be constructed: " << e.what()
                << "\n  It reads heps.tab from $GARFIELD_INSTALL/share/Garfield"
                   "/Data. Run\n    source "
                   "$(dirname $(which garfield 2>/dev/null) 2>/dev/null)"
                   "../share/Garfield/setupGarfield.sh\n  or set "
                   "GARFIELD_INSTALL by hand, or use ionisationModel = "
                   "\"heed\".\n";
      return 1;
    }
    auto& bichsel = *bichselOwner;
    bichsel.SetParticle("pi");
    bichsel.SetMomentum(180.e9);
    if (!bichsel.ComputeCrossSection()) {
      std::cerr << "TrackBichsel::ComputeCrossSection failed.\n";
      return 1;
    }
    // Heed transports Bichsel's delta electrons; NewTrack is never called,
    // so initialise explicitly.
    if (!track.Initialise(&si)) {
      std::cerr << "TrackHeed::Initialise(silicon) failed; cannot transport "
                   "the delta electrons Bichsel produces.\n";
      return 1;
    }
    std::cout << "ionisation: TrackBichsel (silicon cross-sections) + Heed "
                 "delta transport\n  cluster density = "
              << bichsel.GetClusterDensity() << " /cm, dE/dx = "
              << 1.e-3 * bichsel.GetStoppingPower() << " keV/cm\n";
  } else {
    std::cout << "ionisation: TrackHeed on silicon\n";
  }

  // Fills {x, y, z, t} per conduction electron; returns the cluster count.
  auto generatePrimaries =
      [&](const double xStartCm, const double yStartCm,
          std::vector<std::array<double, 4>>& primaries) -> unsigned long {
    primaries.clear();
    unsigned long nClusters = 0;
    if (useBichsel) {
      auto& bichsel = *bichselOwner;
      bichsel.NewTrack(xStartCm, yStartCm, 0., 0., 0., 1., 0.);
      for (const auto& cluster : bichsel.GetClusters()) {
        ++nClusters;
        for (const auto& e : cluster.electrons) {
          const auto cld = track.TransportDeltaElectron(
              e.x, e.y, e.z, e.t, e.e, e.dx, e.dy, e.dz);
          for (const auto& ce : cld.electrons)
            primaries.push_back({ce.x, ce.y, ce.z, ce.t});
        }
      }
    } else {
      track.NewTrack(xStartCm, yStartCm, 0., 0., 0., 1., 0.);
      double xc = 0., yc = 0., zc = 0., tc = 0., ec = 0., extra = 0.;
      int nc = 0;
      while (track.GetCluster(xc, yc, zc, tc, nc, ec, extra)) {
        ++nClusters;
        for (int k = 0; k < nc; ++k) primaries.push_back({xc, yc, zc, tc});
      }
    }
    return nClusters;
  };

  AvalanchePassConfig pc;
  pc.fineStepCm = fineStepCm;
  pc.bulkStepCm = bulkStepCm;
  pc.yFineLoCm = yFineLo;
  pc.yFineHiCm = yFineHi;
  pc.timeWindowNs = cfg.driftWindowNs;
  pc.sizeCap = sizeCap;
  pc.xMin = bx0; pc.yMin = yCollect; pc.zMin = -5.e-4;
  pc.xMax = bx1; pc.yMax = yBot;           pc.zMax = 5.e-4;
  pc.tStart = 0.;
  pc.tStep = cfg.signalStepNs;
  pc.nBins = signalBins;
  pc.enableSignal = true;
  pc.ompDynamicChunk = std::max(1, cfg.avalancheOmpChunk);
  pc.forceSerial = cfg.mipForceSerialAvalanches;
  pc.enableDelayedSignal = dynamicWeightingActive;
  pc.requireDelayedWeightingData = true;
  // Only used for the truncation warning; Sensor reads the component's times.
  pc.delayedSignalTimesNs =
      dynamicWeightingActive ? activeWeightingTimes : std::vector<double>{};
  pc.delayedSignalAveragingOrder = cfg.delayedSignalAveragingOrder;
  pc.silenceGarfieldSetupStdout = cfg.silenceGarfieldSensorSetup;

  ComponentPoisson2d scField;
  SpaceChargeConfig scCfg;
  scCfg.enabled = runScreened;
  scCfg.xMinCm = bx0; scCfg.xMaxCm = bx1;
  scCfg.yMinCm = yTop; scCfg.yMaxCm = yBot;
  scCfg.zMinCm = -5.e-4; scCfg.zMaxCm = 5.e-4;
  scCfg.hMinCm = 0.02e-4;
  scCfg.hMaxCm = 0.50e-4;
  scCfg.zExtentCm = cfg.scZExtentUm * 1.e-4;
  scCfg.tWindowNs = cfg.driftWindowNs;
  scCfg.xProbeCm = x0;
  scCfg.yProbeCm = yGain;
  scCfg.maxIter = cfg.scMaxIter;
  scCfg.relaxation = cfg.scRelaxation;
  scCfg.mixNx = cfg.scMixNx;
  scCfg.mixNy = cfg.scMixNy;
  scCfg.gainTol = cfg.scGainTol;
  scCfg.fieldTol = cfg.scFieldTol;
  scCfg.stableIterations = cfg.scStableIterations;
  scCfg.fieldSampleNx = cfg.scFieldSampleNx;
  scCfg.fieldSampleNy = cfg.scFieldSampleNy;
  scCfg.fieldSampleXHalfWidthCm =
      cfg.scFieldSampleXHalfWidthUm * 1.e-4;
  scCfg.fieldSampleYHalfWidthCm =
      cfg.scFieldSampleYHalfWidthUm * 1.e-4;
  scCfg.verbose = cfg.scConsoleVerbosity >= 2;
  if (scCfg.enabled && !SetupSpaceCharge(scField, scCfg, &si)) {
    std::cerr << "space-charge setup failed; screened MIPs disabled.\n";
    scCfg.enabled = false;
  }
  const bool actuallyRunScreened = runScreened && scCfg.enabled;
  if (cfg.useMeasuredLz && actuallyRunScreened) {
    if (!cfg.measureLz) {
      std::cerr << "useMeasuredLz requires measureLz = true.\n";
      return 3;
    }
    if (!runStatic) {
      std::cerr << "useMeasuredLz takes L_z from the static pass, but "
                   "mipRunStatic is false. Enable the static pass or set "
                   "useMeasuredLz = false to use scZExtentUm.\n";
      return 3;
    }
  }
  if (!runStatic && !actuallyRunScreened) return 2;

  const std::string mipSuffix = biasLabel == "NA" ? "" : "_" + biasLabel + "V";
  const std::string eventsPath = cfg.outDir + "mip_events" + mipSuffix + ".csv";
  const std::string pairedPath = cfg.outDir + "mip_paired" + mipSuffix + ".csv";
  const std::string primariesPath = cfg.outDir + "mip_primaries" + mipSuffix + ".csv";
  const std::string samplesPath = cfg.outDir + "mip_final_samples" + mipSuffix + ".csv";
  const std::string sampleSummaryPath =
      cfg.outDir + "mip_final_sample_summary" + mipSuffix + ".csv";

  std::ofstream fPrimaries;
  if (cfg.mipWritePrimaries) {
    fPrimaries.open(primariesPath);
    fPrimaries << "event,pair,x_cm,y_cm,z_cm,t_ns\n";
  }

  std::ofstream fEvents(eventsPath);
  fEvents << "event,bias,model,mode,nClusters,nPrimary,nTotal,countingGain,"
             "chargeGainPromptAbs,peakOn,peakOff,intOnSigned,intOffSigned,"
             "areaOnPositive,areaOnNegative,areaOnAbsolute,"
             "zeroCrossingOn_ns,lzUsedUm,"
             "scIterations,scConverged,scCharge_e_per_cm,scField_Vcm,"
             "elapsed_s";
  for (const auto& strip : strips) {
    fEvents << "," << strip.label << "_signedAreaOn_arbNs," << strip.label
            << "_signedAreaOff_arbNs," << strip.label
            << "_absoluteAreaFractionOn";
  }
  fEvents << "\n";

  std::ofstream fPaired;
  if (runStatic && actuallyRunScreened) {
    fPaired.open(pairedPath);
    fPaired << "event,bias,model,nClusters,nPrimary,staticCountingGain,"
               "screenedCountingGain,countingRatio,staticChargeGain,"
               "screenedChargeGain,chargeRatio,scIterations,scConverged,"
               "scCharge_e_per_cm,scField_Vcm\n";
  }

  std::ofstream fSamples;
  std::ofstream fSampleSummary;
  if (cfg.mipFinalSampleCount > 0) {
    fSamples.open(samplesPath);
    fSamples << "event,bias,model,mode,sample,nClusters,nPrimary,"
                "nTotal,countingGain,signalEnabled,chargeGainVsCanonicalOff,"
                "scField_Vcm,elapsed_s\n";
    fSampleSummary.open(sampleSummaryPath);
    fSampleSummary << "event,bias,model,mode,nSamples,nPrimary,"
                      "meanCountingGain,sdCountingGain,semCountingGain,"
                      "meanChargeGainVsCanonicalOff,sdChargeGainVsCanonicalOff,"
                      "semChargeGainVsCanonicalOff,scField_Vcm\n";
  }

  std::ofstream iterCsv;
  if (actuallyRunScreened && cfg.scWriteIterationHistory) {
    iterCsv.open(cfg.outDir + "mip_spacecharge_iterations" + mipSuffix
                 + ".csv");
    iterCsv << "event,iteration,relaxation,nPrimary,nTotal,gain,"
               "inputProbeField_Vcm,outputProbeField_Vcm,"
               "freshSignedCharge_e_per_cm,freshAbsCharge_e_per_cm,"
               "mixedSignedCharge_e_per_cm,mixedAbsCharge_e_per_cm,"
               "fieldRelL2,fieldRelMax,fieldSamples,relGainChange,"
               "stableCount,converged\n";
  }

  {
    std::ofstream meta(cfg.outDir + "mip_run_config" + mipSuffix + ".txt");
    meta << "bias=" << biasLabel << "\nmodel=" << cfg.model
         << "\nnMips=" << cfg.nMips
         << "\neventOffset=" << cfg.mipEventOffset
         << "\nrunStatic=" << (runStatic ? 1 : 0)
         << "\nrunScreened=" << (actuallyRunScreened ? 1 : 0)
         << "\nspaceChargeMaster=" << (cfg.spaceCharge ? 1 : 0)
         << "\nscZExtentUm=" << cfg.scZExtentUm
         << "\nuseMeasuredLz=" << (cfg.useMeasuredLz ? 1 : 0)
         << "\nscMaxIter=" << cfg.scMaxIter
         << "\nscRelaxation=" << cfg.scRelaxation
         << "\nscConsoleVerbosity=" << cfg.scConsoleVerbosity
         << "\nscGainTol=" << cfg.scGainTol
         << "\nscConvergenceCriterion=fieldRelL2"
         << "\nscFieldTol=" << cfg.scFieldTol
         << "\nscStableIterations=" << cfg.scStableIterations
         << "\nscMixNx=" << cfg.scMixNx
         << "\nscMixNy=" << cfg.scMixNy
         << "\nscFieldSampleNx=" << cfg.scFieldSampleNx
         << "\nscFieldSampleNy=" << cfg.scFieldSampleNy
         << "\nscFieldSampleXHalfWidthUm="
         << cfg.scFieldSampleXHalfWidthUm
         << "\nscFieldSampleYHalfWidthUm="
         << cfg.scFieldSampleYHalfWidthUm
         << "\nscWriteFinalFieldProfile="
         << (cfg.scWriteFinalFieldProfile ? 1 : 0)
         << "\nscFinalFieldProfilePoints="
         << cfg.scFinalFieldProfilePoints
         << "\nscFinalEvaluationPass="
         << (cfg.scFinalEvaluationPass ? 1 : 0)
         << "\nscWarmStartEvents=" << (cfg.scWarmStartEvents ? 1 : 0)
         << "\nscWarmStartScaleMin=" << cfg.scWarmStartScaleMin
         << "\nscWarmStartScaleMax=" << cfg.scWarmStartScaleMax
         << "\nscPrintStageTimings=" << (cfg.scPrintStageTimings ? 1 : 0)
         << "\navalancheOmpChunk=" << cfg.avalancheOmpChunk
         << "\nionisationModel=" << cfg.ionisationModel
         << "\nrngMode=continuousGarfieldStream"
         << "\nmipForceSerialAvalanches="
         << (cfg.mipForceSerialAvalanches ? 1 : 0)
         << "\nuseResistiveGridWeighting="
         << (cfg.useResistiveGridWeighting ? 1 : 0)
         << "\nmeasureLz=" << (cfg.measureLz ? 1 : 0)
         << "\nlzBandHalfWidthUm=" << cfg.lzBandHalfWidthUm
         << "\nsheetResistanceOhmSq=" << cfg.sheetResistanceOhmSq
         << "\ncouplingEpsR=" << cfg.couplingEpsR
         << "\ncouplingOxideThicknessUm=" << cfg.couplingOxideThicknessUm
         << "\ndcContactUm=" << cfg.dcContactX0Um << ".."
         << cfg.dcContactX1Um
         << "\ngroundedPadsUm=" << [&] {
              std::ostringstream s;
              for (const auto& g : cfg.groundedPadsUm)
                s << g.first << ".." << g.second << ";";
              return s.str();
            }()
         << "\nsheetNodes=" << cfg.sheetNodes
         << "\nsheetTimeSlices=" << cfg.resistiveWeightingTimesNs.size()
         << "\nsheetLastTimeNs=" << cfg.resistiveWeightingTimesNs.back()
         << "\ndepletionEdgeUm=" << cfg.depletionEdgeUm
         << "\nyCollectCm=" << yCollect
         << "\nuseDynamicTcadWeighting="
         << (cfg.useDynamicTcadWeighting ? 1 : 0)
         << "\ndynamicWeightingActive="
         << (dynamicWeightingActive ? 1 : 0)
         << "\nweightingGridFile=" << cfg.weightingGridFile
         << "\nweightingBiasDataFile=" << cfg.weightingBiasDataFile
         << "\nweightingDeltaV=" << cfg.weightingDeltaV
         << "\ndynamicWeightingUsesPotential="
         << (cfg.dynamicWeightingUsesPotential ? 1 : 0)
         << "\nsignalStepNs=" << cfg.signalStepNs
         << "\nsignalWindowNs=" << signalWindowNs
         << "\nsignalBins=" << signalBins
         << "\nenableDelayedSignal="
         << (dynamicWeightingActive ? 1 : 0)
         << "\ndelayedSignalAveragingOrder="
         << cfg.delayedSignalAveragingOrder
         << "\nmipFinalSampleCount=" << cfg.mipFinalSampleCount
         << "\nmipFinalSampleEnableSignal="
         << (cfg.mipFinalSampleEnableSignal ? 1 : 0)
         << "\nmipWritePrimaries=" << (cfg.mipWritePrimaries ? 1 : 0)
         << "\nmipWriteOverlaySignals="
         << (cfg.mipWriteOverlaySignals ? 1 : 0)
         << "\nfineStepNm=" << cfg.fineStepNm
         << "\nbulkStepNm=" << cfg.bulkStepNm
         << "\ndriftWindowNs=" << cfg.driftWindowNs << "\n";
  }

  // one column per MIP, one CSV per strip/pass
  std::vector<int> overlayEventIds;
  std::vector<unsigned long> overlayPrimaryCounts;
  std::vector<std::vector<std::vector<double>>> sigOffOverlay(strips.size());
  std::vector<std::vector<std::vector<double>>> sigStaticOverlay(strips.size());
  std::vector<std::vector<std::vector<double>>> sigScreenedOverlay(strips.size());
  std::vector<double> ensemblePrimaryCounts;
  std::vector<double> ensembleStaticCounting;
  std::vector<double> ensembleStaticCharge;
  std::vector<double> ensembleScreenedCounting;
  std::vector<double> ensembleScreenedCharge;
  std::vector<double> ensembleCountingRatio;
  std::vector<double> ensembleChargeRatio;

  struct ModeResult {
    bool ran = false;
    unsigned long nTotal = 0;
    double countingGain = 0.;
    double chargeGain = 0.;
    SignalMetrics metrics;
    SignalMetrics promptMetrics;
    SignalMetrics delayedMetrics;
    std::vector<std::vector<double>> signal;
    std::vector<std::vector<double>> promptSignal;
    std::vector<std::vector<double>> delayedSignal;
    int scIterations = 0;
    bool scConverged = false;
    ScreeningFieldSample scSample;
    double elapsedS = 0.;
  };

  double lzUsedUm = cfg.scZExtentUm;
  const auto writeModeRow = [&](const int eventId, const char* mode,
                                const unsigned long nClusters,
                                const unsigned long nPrimary,
                                const SignalMetrics& offMetrics,
                                const ModeResult& r) {
    fEvents << eventId << "," << biasLabel << "," << cfg.model << ","
            << mode << "," << nClusters << "," << nPrimary << ","
            << r.nTotal << "," << r.countingGain << "," << r.chargeGain
            << "," << r.metrics.peak[r.metrics.peakStrip] << ","
            << offMetrics.peak[offMetrics.peakStrip] << ","
            << r.metrics.totalIntegral << "," << offMetrics.totalIntegral
            << "," << r.metrics.totalPositiveArea
            << "," << r.metrics.totalNegativeArea
            << "," << r.metrics.totalAbsoluteArea
            << "," << r.metrics.zeroCrossingAfterPeakNs[r.metrics.peakStrip]
            << "," << (std::string(mode) == "screened" ? lzUsedUm : 0.)
            << "," << r.scIterations << "," << (r.scConverged ? 1 : 0)
            << "," << r.scSample.depositedChargeEPerCm << ","
            << r.scSample.magnitudeVcm << "," << r.elapsedS;
    double qAbs = 0.;
    for (const double q : r.metrics.absoluteArea) qAbs += q;
    for (std::size_t k = 0; k < strips.size(); ++k) {
      const double frac = qAbs > 0. ? r.metrics.absoluteArea[k] / qAbs : 0.;
      fEvents << "," << r.metrics.integral[k] << ","
              << offMetrics.integral[k] << "," << frac;
    }
    fEvents << "\n";
  };

  ChargeGrid previousConvergedGrid;
  bool havePreviousConvergedGrid = false;
  unsigned long previousConvergedPrimary = 0;
  double previousConvergedLzUm = 0.;
  if (actuallyRunScreened) previousConvergedGrid.Init(scCfg);

  for (int iMip = 0; iMip < nMipRun; ++iMip) {
    const int eventId = cfg.mipEventOffset + iMip;
    const auto tEventStart = Clock::now();
    std::vector<std::array<double, 4>> primaries;
    const unsigned long nClusters =
        generatePrimaries(x0, yTop + 0.03e-4, primaries);
    const unsigned long nPrimary = primaries.size();
    if (cfg.mipWritePrimaries) {
      for (std::size_t ip = 0; ip < primaries.size(); ++ip) {
        const auto& p = primaries[ip];
        fPrimaries << eventId << "," << ip << "," << p[0] << ","
                   << p[1] << "," << p[2] << "," << p[3] << "\n";
      }
    }
    if (nPrimary == 0) {
      std::cerr << "MIP event " << eventId << " produced no primary pairs; "
                   "skipping.\n";
      continue;
    }

    std::cout << "\n[MIP " << eventId << "] " << nClusters << " clusters, "
              << nPrimary << " primary pairs\n";

    AvalanchePassConfig pcOff = pc;
    pcOff.multiplication = false;
    // prompt-only baseline; skips the costly delayed accumulation
    pcOff.enableDelayedSignal = false;
    std::vector<std::vector<double>> sigOff;
    std::vector<std::vector<double>> sigOffPrompt;
    std::vector<std::vector<double>> sigOffDelayed;
    const auto tGainOff = Clock::now();
    RunAvalanchePass(
cmp, *weightingCmp, perStripWeighting, strips, primaries, pcOff, sigOff, nullptr,
                     cfg.mipPairProgressEvery, "GAIN_OFF",
                     nullptr, nullptr, nullptr, nullptr,
                     &sigOffPrompt, &sigOffDelayed);
    const double gainOffS = ElapsedS(tGainOff);
    if (cfg.scPrintStageTimings) {
      std::cout << "  [timing] gain-off transport=" << gainOffS << " s\n";
    }
    const SignalMetrics offMetrics =
        MeasureSignal(sigOff, pc.tStart, pc.tStep);
    const SignalMetrics offPromptMetrics =
        MeasureSignal(sigOffPrompt, pc.tStart, pc.tStep);
    overlayEventIds.push_back(eventId);
    overlayPrimaryCounts.push_back(nPrimary);
    if (cfg.mipWriteOverlaySignals) {
      for (std::size_t k = 0; k < strips.size(); ++k) {
        sigOffOverlay[k].push_back(sigOff[k]);
      }
    }

    std::vector<double> staticSampleCounting;
    std::vector<double> staticSampleCharge;
    std::vector<double> screenedSampleCounting;
    std::vector<double> screenedSampleCharge;

    const auto writeSampleSummary = [&](const char* mode,
                                        const std::vector<double>& counting,
                                        const std::vector<double>& charge,
                                        const double scFieldVcm) {
      if (!fSampleSummary.is_open() || counting.empty()) return;
      const auto cs = ComputeSampleStats(counting);
      const auto qs = ComputeSampleStats(charge);
      fSampleSummary << eventId << "," << biasLabel << "," << cfg.model
                     << "," << mode << "," << cs.n << "," << nPrimary
                     << "," << cs.mean << "," << cs.sd << "," << cs.sem
                     << "," << qs.mean << "," << qs.sd << "," << qs.sem
                     << "," << scFieldVcm << "\n";
      fSampleSummary.flush();
    };

    const auto runFrozenFieldSamples = [&](const char* mode,
                                            Component* screeningField,
                                            std::vector<double>& counting,
                                            std::vector<double>& charge,
                                            const double scFieldVcm) {
      if (cfg.mipFinalSampleCount <= 0) return;
      AvalanchePassConfig pcSample = pc;
      pcSample.multiplication = true;
      pcSample.enableSignal = cfg.mipFinalSampleEnableSignal;
      for (int is = 0; is < cfg.mipFinalSampleCount; ++is) {
        std::vector<std::vector<double>> sampleSignal;
        std::vector<std::vector<double>> samplePrompt;
        std::vector<std::vector<double>> sampleDelayed;
        const auto ts = Clock::now();
        const unsigned long nTotal = RunAvalanchePass(
cmp, *weightingCmp, perStripWeighting, strips, primaries, pcSample, sampleSignal, nullptr,
            0, std::string(mode) + "_SAMPLE", screeningField,
            nullptr, nullptr, nullptr, &samplePrompt, &sampleDelayed);
        const double gCount = double(nTotal) / nPrimary;
        double gCharge = std::numeric_limits<double>::quiet_NaN();
        if (cfg.mipFinalSampleEnableSignal) {
          const auto promptMetrics =
              MeasureSignal(samplePrompt, pc.tStart, pc.tStep);
          if (offPromptMetrics.totalAbsoluteArea > 0.) {
            gCharge = promptMetrics.totalAbsoluteArea /
                      offPromptMetrics.totalAbsoluteArea;
          }
        }
        counting.push_back(gCount);
        if (std::isfinite(gCharge)) charge.push_back(gCharge);
        if (fSamples.is_open()) {
          fSamples << eventId << "," << biasLabel << "," << cfg.model
                   << "," << mode << "," << is << ","
                   << nClusters << "," << nPrimary << "," << nTotal << ","
                   << gCount << ","
                   << (cfg.mipFinalSampleEnableSignal ? 1 : 0) << ","
                   << gCharge << "," << scFieldVcm << ","
                   << ElapsedS(ts) << "\n";
          fSamples.flush();
        }
      }
      writeSampleSummary(mode, counting, charge, scFieldVcm);
    };

    double measuredLzUm = 0.;  // 0 = no valid measurement
    double measuredSigmaZUm = 0.;
    ModeResult staticResult;
    if (runStatic) {
      const auto t0 = Clock::now();
      staticResult.ran = true;
      AvalanchePassConfig pcStatic = pc;
      pcStatic.multiplication = true;
      std::ofstream pairs;
      std::ofstream* pairPtr = nullptr;
      if (cfg.mipWritePairFiles) {
        const std::string path = cfg.outDir + "mip_pairs" + mipSuffix + "_"
            + EventTag(eventId) + "_static.csv";
        pairs.open(path);
        pairs << "x_um,y_um,ne\n";
        pairPtr = &pairs;
      }
      std::ofstream driftFs;
      if (cfg.mipWriteDriftLines && iMip == 0) {
        driftFs.open(cfg.outDir + "driftlines" + mipSuffix + ".csv");
        driftFs << "carrier,species,x_um,y_um,t_ns\n";
        pcStatic.driftCsv = &driftFs;
        pcStatic.driftMaxCarriers = cfg.driftMaxCarriers;
      }
      TransverseExtentAccumulator lzGain;
      staticResult.nTotal = RunAvalanchePass(
cmp, *weightingCmp, perStripWeighting, strips, primaries, pcStatic, staticResult.signal, pairPtr,
          cfg.mipPairProgressEvery, "STATIC",
          nullptr, nullptr, nullptr, nullptr,
          &staticResult.promptSignal, &staticResult.delayedSignal,
          cfg.measureLz ? &lzGain : nullptr,
          yGain - cfg.lzBandHalfWidthUm * 1.e-4,
          yGain + cfg.lzBandHalfWidthUm * 1.e-4);
      if (cfg.measureLz) {
        const auto teGain = lzGain.Finish();
        if (teGain.valid) {
          measuredLzUm = teGain.lzCm * 1.e4;
          measuredSigmaZUm = teGain.sigmaZCm * 1.e4;
          std::cout << "  [L_z] gain-layer band +/-" << cfg.lzBandHalfWidthUm
                    << " um about y=" << yGain * 1.e4 << " um: sigma_z = "
                    << measuredSigmaZUm << " um,  L_z = sqrt(12)*sigma_z = "
                    << measuredLzUm << " um   (" << teGain.nSegments
                    << " segments, " << teGain.weightNsE << " carrier-ns)\n"
                    << "  [L_z] legacy fitted scZExtentUm = "
                    << cfg.scZExtentUm << " um  -> measured/fitted = "
                    << (measuredLzUm / cfg.scZExtentUm) << "\n";
        } else {
          measuredLzUm = 0.;
          std::cerr << "  [L_z] measurement invalid ("
                    << teGain.nCarriersNoPath
                    << " carriers had no stored path)\n";
        }
      }
      staticResult.metrics =
          MeasureSignal(staticResult.signal, pc.tStart, pc.tStep);
      staticResult.promptMetrics =
          MeasureSignal(staticResult.promptSignal, pc.tStart, pc.tStep);
      staticResult.delayedMetrics =
          MeasureSignal(staticResult.delayedSignal, pc.tStart, pc.tStep);
      staticResult.countingGain = double(staticResult.nTotal) / nPrimary;
      // prompt |area|: AC-coupled signed lobes cancel
      staticResult.chargeGain = offPromptMetrics.totalAbsoluteArea > 0.
          ? staticResult.promptMetrics.totalAbsoluteArea /
            offPromptMetrics.totalAbsoluteArea : 0.;
      staticResult.elapsedS = ElapsedS(t0);
      if (cfg.mipWriteOverlaySignals) {
        for (std::size_t k = 0; k < strips.size(); ++k) {
          sigStaticOverlay[k].push_back(staticResult.signal[k]);
        }
      }
      writeModeRow(eventId, "static", nClusters, nPrimary, offMetrics,
                   staticResult);
      std::cout << "  static: Gcount=" << staticResult.countingGain
                << " Gcharge=" << staticResult.chargeGain << "\n";
      runFrozenFieldSamples("static", nullptr,
                            staticSampleCounting, staticSampleCharge, 0.);
    }

    lzUsedUm = cfg.scZExtentUm;
    if (cfg.useMeasuredLz && actuallyRunScreened) {
      if (measuredLzUm > 0.) {
        lzUsedUm = measuredLzUm;
      } else {
        std::cerr << "  [L_z] no valid measurement this event; falling back "
                     "to scZExtentUm = " << cfg.scZExtentUm << " um\n";
      }
    }
    scCfg.zExtentCm = lzUsedUm * 1.e-4;
    if (actuallyRunScreened) {
      std::cout << "  [L_z] screening will use L_z = " << lzUsedUm << " um"
                << (cfg.useMeasuredLz && measuredLzUm > 0. ? " (measured)"
                                                           : " (configured)")
                << "\n";
    }

    ModeResult screenedResult;
    if (actuallyRunScreened) {
      const auto t0 = Clock::now();
      screenedResult.ran = true;

      double gainPrev = std::numeric_limits<double>::quiet_NaN();
      ScreeningFieldGridSample fieldGridPrev;
      bool haveFieldGridPrev = false;
      int stableCount = 0;
      std::vector<std::vector<double>> sigIter;
      unsigned long nIterTotal = 0;

      ChargeGrid mixed, fresh;
      mixed.Init(scCfg);
      fresh.Init(scCfg);

      bool usedWarmStart = false;
      if (cfg.scWarmStartEvents && havePreviousConvergedGrid &&
          previousConvergedPrimary > 0) {
        mixed = previousConvergedGrid;
        // rho is per unit L_z, so rescale by the L_z ratio too
        const double lzScale = (previousConvergedLzUm > 0. && lzUsedUm > 0.)
            ? previousConvergedLzUm / lzUsedUm : 1.;
        const double rawScale = lzScale * static_cast<double>(nPrimary) /
                                previousConvergedPrimary;
        const double scaleLo = std::min(cfg.scWarmStartScaleMin,
                                        cfg.scWarmStartScaleMax);
        const double scaleHi = std::max(cfg.scWarmStartScaleMin,
                                        cfg.scWarmStartScaleMax);
        const double scale = std::clamp(rawScale, scaleLo, scaleHi);
        for (double& q : mixed.rho) q *= scale;
        LoadChargeGrid(scField, mixed);
        if (scField.Solve()) {
          usedWarmStart = true;
          fieldGridPrev = SampleScreeningFieldGrid(scField, scCfg);
          haveFieldGridPrev = true;
          if (cfg.scConsoleVerbosity >= 1) {
            std::cout << "  [spacecharge] warm start from previous event"
                      << " (primary scale=" << scale << ")\n";
          }
        } else {
          std::cerr << "  [spacecharge] warm-start solve failed; "
                       "falling back to zero.\n";
          mixed.Clear();
        }
      }
      if (!usedWarmStart) {
        scField.ClearCharge();
        if (!scField.Solve()) {
          std::cerr << "  [spacecharge] zero-field reset failed.\n";
        }
      }

      AvalanchePassConfig pcIter = pc;
      pcIter.enableSignal = false;  // intermediate signals are discarded
      bool haveConvergedIterationSignal = false;
      unsigned long convergedIterationTotal = 0;
      std::vector<std::vector<double>> convergedIterationSignal;

      for (int it = 0; it < scCfg.maxIter; ++it) {
        if (cfg.scConsoleVerbosity >= 2) {
          std::cout << "  [spacecharge] iteration " << it + 1 << "/"
                    << scCfg.maxIter << "\n";
        }
        const auto inputFieldSample = SampleScreeningField(scField, scCfg);
        fresh.Clear();
        // Without a final pass, record the signal on the likely last iteration.
        pcIter.enableSignal = !cfg.scFinalEvaluationPass &&
            !cfg.mipWritePairFiles &&
            stableCount >= std::max(0, scCfg.stableIterations - 1);

        const auto tTransport = Clock::now();
        nIterTotal = RunAvalanchePass(
cmp, *weightingCmp, perStripWeighting, strips, primaries, pcIter, sigIter, nullptr,
            cfg.mipPairProgressEvery, "SCREEN_ITER", &scField, nullptr,
            &scCfg, &fresh);
        const double transportS = ElapsedS(tTransport);

        const auto tApply = Clock::now();
        ApplyRelaxedCharge(scField, mixed, fresh, scCfg.relaxation);
        const double applyS = ElapsedS(tApply);

        const auto tSolve = Clock::now();
        const bool solved = scField.Solve();
        const double solveS = ElapsedS(tSolve);
        if (!solved) {
          std::cerr << "  [spacecharge] Solve() failed; stopping iterations.\n";
          break;
        }

        const auto tSample = Clock::now();
        screenedResult.scIterations = it + 1;
        screenedResult.scSample = SampleScreeningField(scField, scCfg);
        const auto fieldGrid = SampleScreeningFieldGrid(scField, scCfg);
        const double sampleS = ElapsedS(tSample);
        const auto fieldChange = haveFieldGridPrev
            ? CompareScreeningFieldGrids(fieldGridPrev, fieldGrid)
            : ScreeningFieldChange{};
        if (cfg.scConsoleVerbosity >= 2) {
          ReportScreeningField(scField, scCfg);
        }
        const double gain = double(nIterTotal) / nPrimary;
        const double fld = screenedResult.scSample.magnitudeVcm;
        const double dG = std::isfinite(gainPrev)
            ? std::abs(gain - gainPrev) /
              std::max(1.e-12, std::abs(gainPrev))
            : std::numeric_limits<double>::quiet_NaN();
        const double dE = fieldChange.relativeL2;

        // converge on the field; dG is too noisy to gate
        const bool stableNow = std::isfinite(dE) && dE < scCfg.fieldTol;
        stableCount = stableNow ? stableCount + 1 : 0;

        if (cfg.scConsoleVerbosity >= 1) {
          std::cout << "  [SC " << it + 1 << "/" << scCfg.maxIter
                    << "] G=" << gain << " Esc=" << fld << " V/cm"
                    << " dG=" << dG << " dE=" << dE
                    << " stable=" << stableCount << "/"
                    << scCfg.stableIterations << "\n";
        }
        const bool convergedNow =
            stableCount >= std::max(1, scCfg.stableIterations);
        if (cfg.scPrintStageTimings) {
          std::cout << "    [timing] transport=" << transportS
                    << " s apply=" << applyS
                    << " s solve=" << solveS
                    << " s sample=" << sampleS << " s\n";
        }
        if (convergedNow && pcIter.enableSignal) {
          haveConvergedIterationSignal = true;
          convergedIterationTotal = nIterTotal;
          convergedIterationSignal = sigIter;
        }
        if (iterCsv.is_open()) {
          iterCsv << eventId << "," << it + 1 << ","
                  << scCfg.relaxation << "," << nPrimary << ","
                  << nIterTotal << "," << gain << ","
                  << inputFieldSample.magnitudeVcm << "," << fld << ","
                  << fresh.Total() << "," << fresh.AbsoluteTotal() << ","
                  << mixed.Total() << "," << mixed.AbsoluteTotal() << ","
                  << fieldChange.relativeL2 << ","
                  << fieldChange.relativeMax << ","
                  << fieldChange.nCompared << "," << dG << ","
                  << stableCount << "," << (convergedNow ? 1 : 0) << "\n";
          iterCsv.flush();
        }
        if (convergedNow) {
          screenedResult.scConverged = true;
          break;
        }
        gainPrev = gain;
        fieldGridPrev = fieldGrid;
        haveFieldGridPrev = true;
      }

      if (cfg.scWarmStartEvents && screenedResult.scConverged) {
        previousConvergedGrid = mixed;
        previousConvergedPrimary = nPrimary;
        previousConvergedLzUm = lzUsedUm;
        havePreviousConvergedGrid = true;
      }

      std::ofstream pairs;
      std::ofstream* pairPtr = nullptr;
      if (cfg.mipWritePairFiles) {
        const std::string path = cfg.outDir + "mip_pairs" + mipSuffix + "_"
            + EventTag(eventId) + "_screened.csv";
        pairs.open(path);
        pairs << "x_um,y_um,ne\n";
        pairPtr = &pairs;
      }

      AvalanchePassConfig pcFinal = pc;
      pcFinal.multiplication = true;
      pcFinal.enableSignal = true;
      if (!dynamicWeightingActive &&
          !cfg.scFinalEvaluationPass && !cfg.mipWritePairFiles &&
          haveConvergedIterationSignal) {
        screenedResult.nTotal = convergedIterationTotal;
        screenedResult.signal = std::move(convergedIterationSignal);
        screenedResult.promptSignal = screenedResult.signal;
        screenedResult.delayedSignal.assign(
            strips.size(), std::vector<double>(pc.nBins, 0.));
        if (cfg.scPrintStageTimings) {
          std::cout << "  [timing] reused converged-iteration signal; "
                       "skipped SCREEN_FINAL\n";
        }
      } else {
        if (!cfg.scFinalEvaluationPass && cfg.mipWritePairFiles) {
          std::cout << "  [spacecharge] pair-file output requires "
                       "SCREEN_FINAL.\n";
        } else if (!cfg.scFinalEvaluationPass &&
                   !haveConvergedIterationSignal) {
          std::cerr << "  [spacecharge] no converged iteration signal was "
                       "available; running SCREEN_FINAL.\n";
        }
        const auto tFinal = Clock::now();
        screenedResult.nTotal = RunAvalanchePass(
cmp, *weightingCmp, perStripWeighting, strips, primaries, pcFinal, screenedResult.signal,
            pairPtr, cfg.mipPairProgressEvery, "SCREEN_FINAL",
            &scField, nullptr, nullptr, nullptr,
            &screenedResult.promptSignal, &screenedResult.delayedSignal);
        if (cfg.scPrintStageTimings) {
          std::cout << "  [timing] final screened transport="
                    << ElapsedS(tFinal) << " s\n";
        }
      }
      screenedResult.scSample = SampleScreeningField(scField, scCfg);
      if (cfg.scWriteFinalFieldProfile) {
        const std::string fieldPath = cfg.outDir + "spacecharge_field" +
            mipSuffix + "_" + EventTag(eventId) + ".csv";
        DumpCombinedFieldProfile(cmp, scField, x0, yTop, yBot,
                                 cfg.scFinalFieldProfilePoints, fieldPath,
                                 cfg.scConsoleVerbosity >= 1);
      }
      screenedResult.metrics =
          MeasureSignal(screenedResult.signal, pc.tStart, pc.tStep);
      screenedResult.promptMetrics =
          MeasureSignal(screenedResult.promptSignal, pc.tStart, pc.tStep);
      screenedResult.delayedMetrics =
          MeasureSignal(screenedResult.delayedSignal, pc.tStart, pc.tStep);
      screenedResult.countingGain = double(screenedResult.nTotal) / nPrimary;
      screenedResult.chargeGain = offPromptMetrics.totalAbsoluteArea > 0.
          ? screenedResult.promptMetrics.totalAbsoluteArea /
            offPromptMetrics.totalAbsoluteArea : 0.;
      screenedResult.elapsedS = ElapsedS(t0);
      if (cfg.mipWriteOverlaySignals) {
        for (std::size_t k = 0; k < strips.size(); ++k) {
          sigScreenedOverlay[k].push_back(screenedResult.signal[k]);
        }
      }
      writeModeRow(eventId, "screened", nClusters, nPrimary, offMetrics,
                   screenedResult);
      std::cout << "  screened: Gcount=" << screenedResult.countingGain
                << " Gcharge=" << screenedResult.chargeGain
                << " field=" << screenedResult.scSample.magnitudeVcm
                << " V/cm\n";
      runFrozenFieldSamples("screened", &scField,
                            screenedSampleCounting, screenedSampleCharge,
                            screenedResult.scSample.magnitudeVcm);
    }

    ensemblePrimaryCounts.push_back(static_cast<double>(nPrimary));
    if (staticResult.ran) {
      ensembleStaticCounting.push_back(staticResult.countingGain);
      ensembleStaticCharge.push_back(staticResult.chargeGain);
    }
    if (screenedResult.ran) {
      ensembleScreenedCounting.push_back(screenedResult.countingGain);
      ensembleScreenedCharge.push_back(screenedResult.chargeGain);
    }

    if (staticResult.ran && screenedResult.ran) {
      const double countRatio = staticResult.countingGain != 0.
          ? screenedResult.countingGain / staticResult.countingGain : 0.;
      const double chargeRatio = staticResult.chargeGain != 0.
          ? screenedResult.chargeGain / staticResult.chargeGain : 0.;
      ensembleCountingRatio.push_back(countRatio);
      ensembleChargeRatio.push_back(chargeRatio);
      fPaired << eventId << "," << biasLabel << "," << cfg.model << ","
              << nClusters << "," << nPrimary << ","
              << staticResult.countingGain << ","
              << screenedResult.countingGain << "," << countRatio << ","
              << staticResult.chargeGain << ","
              << screenedResult.chargeGain << "," << chargeRatio << ","
              << screenedResult.scIterations << ","
              << (screenedResult.scConverged ? 1 : 0) << ","
              << screenedResult.scSample.depositedChargeEPerCm << ","
              << screenedResult.scSample.magnitudeVcm << "\n";
    }

    if (cfg.mipWritePerEventSignals) {
      std::ofstream fs(cfg.outDir + "signal" + mipSuffix + "_"
                       + EventTag(eventId) + ".csv");
      fs << "bin,t_ns";
      for (const auto& strip : strips) {
        fs << "," << strip.label << "_GAIN_OFF_TOTAL"
           << "," << strip.label << "_GAIN_OFF_PROMPT"
           << "," << strip.label << "_GAIN_OFF_DELAYED";
        if (staticResult.ran) {
          fs << "," << strip.label << "_STATIC_TOTAL"
             << "," << strip.label << "_STATIC_PROMPT"
             << "," << strip.label << "_STATIC_DELAYED";
        }
        if (screenedResult.ran) {
          fs << "," << strip.label << "_SCREENED_TOTAL"
             << "," << strip.label << "_SCREENED_PROMPT"
             << "," << strip.label << "_SCREENED_DELAYED";
        }
      }
      fs << "\n";
      for (unsigned int b = 0; b < pc.nBins; ++b) {
        fs << b << "," << (b + 0.5) * pc.tStep;
        for (std::size_t k = 0; k < strips.size(); ++k) {
          fs << "," << sigOff[k][b]
             << "," << sigOffPrompt[k][b]
             << "," << sigOffDelayed[k][b];
          if (staticResult.ran) {
            fs << "," << staticResult.signal[k][b]
               << "," << staticResult.promptSignal[k][b]
               << "," << staticResult.delayedSignal[k][b];
          }
          if (screenedResult.ran) {
            fs << "," << screenedResult.signal[k][b]
               << "," << screenedResult.promptSignal[k][b]
               << "," << screenedResult.delayedSignal[k][b];
          }
        }
        fs << "\n";
      }
    }

    if (cfg.mipProgressEvery > 0 &&
        (iMip + 1) % cfg.mipProgressEvery == 0) {
      std::cout << "[MIP ensemble] completed " << iMip + 1 << "/"
                << cfg.nMips << " events in " << ElapsedS(tEventStart)
                << " s for this event\n";
    }
  }

  if (cfg.mipWriteOverlaySignals && !overlayEventIds.empty()) {
    const auto writeOverlay = [&](const std::string& mode,
                                  const std::vector<std::vector<std::vector<double>>>& data) {
      for (std::size_t k = 0; k < strips.size(); ++k) {
        if (data[k].size() != overlayEventIds.size()) {
          std::cerr << "overlay " << mode << "/" << strips[k].label
                    << " has " << data[k].size() << " signals for "
                    << overlayEventIds.size() << " events; skipping file.\n";
          continue;
        }
        const std::string path = cfg.outDir + "signal_overlay_" + mode + "_"
            + strips[k].label + mipSuffix + ".csv";
        std::ofstream out(path);
        out << "time_ns";
        for (const int eventId : overlayEventIds) out << ",trk" << eventId;
        out << "\n";
        out << "nPairs";
        for (const auto nPrimary : overlayPrimaryCounts) out << "," << nPrimary;
        out << "\n";
        for (unsigned int b = 0; b < pc.nBins; ++b) {
          out << (b + 0.5) * pc.tStep;
          for (std::size_t i = 0; i < overlayEventIds.size(); ++i) {
            out << "," << data[k][i][b];
          }
          out << "\n";
        }
        std::cout << "wrote " << path << "\n";
      }
    };

    writeOverlay("gainOff", sigOffOverlay);
    if (runStatic) writeOverlay("static", sigStaticOverlay);
    if (actuallyRunScreened) writeOverlay("screened", sigScreenedOverlay);
  }

  const std::string ensembleSummaryPath =
      cfg.outDir + "mip_ensemble_summary" + mipSuffix + ".csv";
  {
    std::ofstream out(ensembleSummaryPath);
    out << "quantity,n,mean,sd,sem\n";
    const auto writeStats = [&](const char* name,
                                const std::vector<double>& values) {
      if (values.empty()) return;
      const auto st = ComputeSampleStats(values);
      out << name << "," << st.n << "," << st.mean << ","
          << st.sd << "," << st.sem << "\n";
    };
    writeStats("nPrimary", ensemblePrimaryCounts);
    writeStats("staticCountingGain", ensembleStaticCounting);
    writeStats("staticPromptAbsGain", ensembleStaticCharge);
    writeStats("screenedCountingGain", ensembleScreenedCounting);
    writeStats("screenedPromptAbsGain", ensembleScreenedCharge);
    writeStats("screenedOverStaticCounting", ensembleCountingRatio);
    writeStats("screenedOverStaticPromptAbs", ensembleChargeRatio);
  }

  std::cout << "wrote " << ensembleSummaryPath << "\n";
  std::cout << "wrote " << eventsPath << "\n";
  if (cfg.mipWritePrimaries) std::cout << "wrote " << primariesPath << "\n";
  if (runStatic && actuallyRunScreened) std::cout << "wrote " << pairedPath << "\n";
  std::cout << "[timer] total run so far: " << ElapsedS(tRunStart)
            << " s\n";

  if (cfg.doPositionScan) {
    const std::string scanPath =
        cfg.outDir + "position_scan" + mipSuffix + ".csv";
    std::ofstream fScan(scanPath);
    fScan << "x_track_um,event,nPrimary,gain";
    for (const auto& s2 : strips) fScan << "," << s2.label << "_peak_uA";
    for (const auto& s2 : strips) fScan << "," << s2.label << "_charge_fC";
    for (const auto& s2 : strips) fScan << "," << s2.label << "_peakNorm";
    for (const auto& s2 : strips) fScan << "," << s2.label << "_chargeNorm";
    for (const auto& s2 : strips) fScan << "," << s2.label << "_tPeak_ns";
    for (const auto& s2 : strips) fScan << "," << s2.label << "_tCfd_ns";
    for (const auto& s2 : strips) fScan << "," << s2.label << "_qPrompt_fC";
    fScan << ",x_reco_um,dx_um,tSeed_ns,seedStrip\n";

    const int nPts = cfg.scanRandomUniform
        ? cfg.scanNRandomEvents
        : static_cast<int>(std::floor((cfg.scanXEndUm - cfg.scanXStartUm) /
                                      cfg.scanXStepUm)) + 1;
    const int nPer = cfg.scanRandomUniform ? 1 : cfg.scanMipsPerPoint;
    // separate stream so scan size doesn't shift the MIP sequence
    std::mt19937_64 scanRng(cfg.scanRandomSeed);
    std::uniform_real_distribution<double> scanX(cfg.scanXStartUm,
                                                 cfg.scanXEndUm);
    std::cout << "\n=== position scan: "
              << (cfg.scanRandomUniform ? "UNIFORM RANDOM " : "grid ")
              << nPts << (cfg.scanRandomUniform ? " throws" : " points")
              << " x " << nPer << " events, gain "
              << (cfg.scanUseGain ? "ON" : "OFF") << " ===" << std::endl;
    const auto tScan = Clock::now();
    int nRun = 0;

    for (int ip = 0; ip < nPts; ++ip) {
      const double xs = cfg.scanRandomUniform
          ? scanX(scanRng)
          : cfg.scanXStartUm + ip * cfg.scanXStepUm;
      // 0.13 um offset keeps the track off mesh lines
      const double xTrackCm = xs * 1.e-4 + 0.13e-4;
      if (xTrackCm <= bx0 || xTrackCm >= bx1) {
        std::cout << "  x = " << xs << " um lies outside the map; skipped\n";
        continue;
      }
      for (int ie = 0; ie < nPer; ++ie) {
        std::vector<std::array<double, 4>> prim;
        generatePrimaries(xTrackCm, yTop + 0.03e-4, prim);
        const unsigned long nPrim = prim.size();
        if (nPrim == 0) continue;

        AvalanchePassConfig pcScan = pc;
        pcScan.multiplication = cfg.scanUseGain;
        pcScan.enableSignal = true;
        if (cfg.scanForceSerial) pcScan.forceSerial = true;
        std::vector<std::vector<double>> sig;
        const unsigned long nTot = RunAvalanchePass(
            cmp, *weightingCmp, perStripWeighting, strips, prim, pcScan, sig,
            nullptr, 0, "SCAN");
        const double gain = cfg.scanUseGain
            ? static_cast<double>(nTot) / static_cast<double>(nPrim) : 1.;
        const double norm = static_cast<double>(nPrim) * gain;

        std::vector<double> peak(strips.size(), 0.), chg(strips.size(), 0.);
        std::vector<double> tPeak(strips.size(),
                                  std::numeric_limits<double>::quiet_NaN());
        std::vector<double> tCfd(strips.size(),
                                 std::numeric_limits<double>::quiet_NaN());
        for (std::size_t k = 0; k < strips.size(); ++k) {
          std::size_t iPeak = 0;
          for (std::size_t b = 0; b < sig[k].size(); ++b) {
            const double v = sig[k][b];
            if (std::abs(v) > std::abs(peak[k])) { peak[k] = v; iPeak = b; }
            chg[k] += v * pc.tStep;
          }
          if (peak[k] == 0.) continue;
          tPeak[k] = pc.tStart + (iPeak + 0.5) * pc.tStep;
          // CFD: walk back from the peak (forward scan hits crosstalk lobes)
          const double thr = cfg.scanCfdFraction * peak[k];
          for (std::size_t b = iPeak; b > 0; --b) {
            const double y1 = sig[k][b], y0 = sig[k][b - 1];
            if ((peak[k] > 0. && y0 < thr && y1 >= thr) ||
                (peak[k] < 0. && y0 > thr && y1 <= thr)) {
              const double f = (y1 == y0) ? 0. : (thr - y0) / (y1 - y0);
              tCfd[k] = pc.tStart + (double(b) - 0.5 + f) * pc.tStep;
              break;
            }
          }
        }
        // Integrate to the first zero crossing after the peak; the full
        // AC-coupled integral cancels.
        std::vector<double> qPrompt(strips.size(), 0.);
        for (std::size_t k = 0; k < strips.size(); ++k) {
          if (peak[k] == 0.) continue;
          const double sgn = peak[k] > 0. ? 1. : -1.;
          for (std::size_t b = 0; b < sig[k].size(); ++b) {
            if (sgn * sig[k][b] <= 0. && pc.tStart + b * pc.tStep >
                (std::isnan(tPeak[k]) ? 0. : tPeak[k])) break;
            qPrompt[k] += sig[k][b] * pc.tStep;
          }
        }
        // centroid weighted by peak amplitude
        double sw = 0., swx = 0.;
        std::size_t seed = 0;
        for (std::size_t k = 0; k < strips.size(); ++k) {
          const double w = std::abs(peak[k]);
          sw += w;
          swx += w * strips[k].centerUm;
          if (w > std::abs(peak[seed])) seed = k;
        }
        const double xReco = sw > 0.
            ? swx / sw : std::numeric_limits<double>::quiet_NaN();

        fScan << xs << "," << ie << "," << nPrim << "," << gain;
        for (double v : peak) fScan << "," << v;
        for (double v : chg) fScan << "," << v;
        for (double v : peak) fScan << "," << v / norm;
        for (double v : chg) fScan << "," << v / norm;
        for (double v : tPeak) fScan << "," << v;
        for (double v : tCfd) fScan << "," << v;
        for (double v : qPrompt) fScan << "," << v;
        // Falls back to tPeak: without an amplifier model the seed strip
        // rises as a step, so CFD is usually undefined.
        const double tSeed = std::isnan(tCfd[seed]) ? tPeak[seed] : tCfd[seed];
        fScan << "," << xReco << "," << (xReco - xs) << "," << tSeed
              << "," << seed << "\n";
        fScan.flush();
      }
      ++nRun;
      std::cout << "  x = " << xs << " um done (" << nRun << " points, "
                << ElapsedS(tScan) << " s elapsed)" << std::endl;
    }
    std::cout << "wrote " << scanPath << "  [" << nRun << " points, "
              << ElapsedS(tScan) << " s]" << std::endl;
  }

  std::cout << "[timer] TOTAL RUNTIME: " << ElapsedS(tRunStart)
            << " s" << std::endl;
  std::cout << "Done.\n";
}