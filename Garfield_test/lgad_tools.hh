#pragma once

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <fstream>
#include <limits>
#include <streambuf>
#include <string>
#include <utility>
#include <vector>

#include "Garfield/AvalancheMC.hh"
#include "Garfield/Component.hh"
#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/ComponentPoisson2d.hh"
#include "Garfield/MediumSilicon.hh"

using Clock = std::chrono::steady_clock;

/* Scales the field seen by impact ionisation only (drift is validated, left
   untouched). f = 1 - kappa*G; kappa = 0.00313 for lgad_deepjte/Okuto. */
class ScreenedSilicon : public Garfield::MediumSilicon {
 public:
  void SetFieldScale(const double f) { m_fscale = f; }
  double GetFieldScale() const { return m_fscale; }

  bool ElectronTownsend(const double ex, const double ey, const double ez,
                        const double bx, const double by, const double bz,
                        double& alpha) override {
    return Garfield::MediumSilicon::ElectronTownsend(
        ex * m_fscale, ey * m_fscale, ez * m_fscale, bx, by, bz, alpha);
  }

  bool HoleTownsend(const double ex, const double ey, const double ez,
                    const double bx, const double by, const double bz,
                    double& alpha) override {
    return Garfield::MediumSilicon::HoleTownsend(
        ex * m_fscale, ey * m_fscale, ez * m_fscale, bx, by, bz, alpha);
  }

 private:
  double m_fscale = 1.0;
};


double ElapsedS(const Clock::time_point& t0);

/* Space-charge feedback via ComponentPoisson2d, layered on the TCAD map:
     sensor.AddComponent(&cmp);  // TCAD first: sensor takes its medium from it
     sensor.AddComponent(&sc);   // screening perturbation
   No SetDoping: the TCAD map already holds the depletion charge. */
struct SpaceChargeConfig {
  // FEM region [cm]; edges must be axis-parallel
  double xMinCm = 0., xMaxCm = 0.;
  double yMinCm = 0., yMaxCm = 0.;
  double zMinCm = -5.e-4, zMaxCm = 5.e-4;

  double hMinCm = 0.02e-4;  // element size near electrodes
  double hMaxCm = 0.50e-4;  // element size in bulk

  // Effective z-extent of the charge column [cm]. Physics input: screening
  // scales linearly with 1/zExtentCm.
  double zExtentCm = 1.e-4;

  double tWindowNs = 1.;  // charge time-averaging window

  int maxIter = 6;

  // rho = (1-lambda)*rho_prev + lambda*rho_new; optimum ~0.67, 0.5 for margin
  double relaxation = 0.5;

  // Mixing grid nodes; blend density, not amplitude (cloud shape changes)
  int mixNx = 101, mixNy = 501;

  // Converged when the sampled screening field changes < fieldTol for
  // stableIterations consecutive iterations.
  double fieldTol = 0.01;
  int stableIterations = 3;
  int fieldSampleNx = 11;
  int fieldSampleNy = 41;
  double fieldSampleXHalfWidthCm = 5.e-4;
  double fieldSampleYHalfWidthCm = 1.e-4;

  double gainTol = 0.01;  // diagnostic only

  double xProbeCm = 0., yProbeCm = 0.;

  bool enabled = false;
  bool verbose = true;
};

// Builds the FEM region, grounded electrodes and mesh (no doping).
bool SetupSpaceCharge(Garfield::ComponentPoisson2d& sc,
                      const SpaceChargeConfig& cfg,
                      Garfield::Medium* medium);

/* Measured out-of-plane extent L_z (replaces the fitted 0.21 um).
   Residence-time weighted (w*dt) z spread of transported carriers;
   lz = sqrt(12)*sigma_z. Restrict yLo/yHi to the gain layer. */
struct TransverseExtent {
  double zBarCm = 0.;
  double sigmaZCm = 0.;
  double lzCm = 0.;
  double weightNsE = 0.;  // sum of w*dt [carrier-ns]
  std::size_t nSegments = 0;
  std::size_t nCarriersNoPath = 0;
  bool valid = false;
};

// Requires aval.EnableDriftLines(true). yLoCm >= yHiCm uses all depths.
TransverseExtent MeasureTransverseExtent(const Garfield::AvalancheMC& aval,
                                         double yLoCm = 0., double yHiCm = 0.,
                                         bool useElectrons = true,
                                         bool useHoles = true);

struct TransverseExtentAccumulator {
  double sumW = 0., sumWz = 0., sumWz2 = 0.;
  std::size_t nSegments = 0, nCarriersNoPath = 0;
  void Add(const Garfield::AvalancheMC& aval, double yLoCm = 0.,
           double yHiCm = 0., bool useElectrons = true, bool useHoles = true);
  void Merge(const TransverseExtentAccumulator& o);
  TransverseExtent Finish() const;
};

/* Deposits carriers (e-, h+) weighted by residence time; no Solve().
   Requires aval.EnableDriftLines(true), else silently deposits nothing.
   Returns segments outside the mesh. */
std::size_t DepositAvalancheCharge(const Garfield::AvalancheMC& aval,
                                   Garfield::ComponentPoisson2d& sc,
                                   const SpaceChargeConfig& cfg);

// Cloud-in-cell charge density grid for blending iterations.
struct ChargeGrid {
  int nx = 0, ny = 0;
  double xMinCm = 0., xMaxCm = 0., yMinCm = 0., yMaxCm = 0.;
  double invDx = 0., invDy = 0.;
  std::vector<double> rho;  // [e / cm of depth]

  void Init(const SpaceChargeConfig& cfg);
  void Clear() { std::fill(rho.begin(), rho.end(), 0.); }
  bool Contains(double xCm, double yCm) const noexcept {
    return xCm >= xMinCm && xCm <= xMaxCm &&
           yCm >= yMinCm && yCm <= yMaxCm;
  }
  void Deposit(double xCm, double yCm, double q);
  void DepositUnchecked(double xCm, double yCm, double q) noexcept;
  double Total() const;
  double AbsoluteTotal() const;
  double NodeX(int ix) const;
  double NodeY(int iy) const;
};

// As DepositAvalancheCharge, into the grid. Returns segments outside it.
std::size_t DepositAvalancheIntoGrid(const Garfield::AvalancheMC& aval,
                                     ChargeGrid& grid,
                                     const SpaceChargeConfig& cfg);

// Neither calls Solve().
void LoadChargeGrid(Garfield::ComponentPoisson2d& sc,
                    const ChargeGrid& grid);

void ApplyRelaxedCharge(Garfield::ComponentPoisson2d& sc, ChargeGrid& mixed,
                        const ChargeGrid& fresh, double lambda);

struct ScreeningFieldSample {
  double exVcm = 0.;
  double eyVcm = 0.;
  double magnitudeVcm = 0.;
  double potentialV = 0.;
  double depositedChargeEPerCm = 0.;
  int status = 0;
};

struct ScreeningFieldGridSample {
  std::vector<double> exVcm;
  std::vector<double> eyVcm;
  std::size_t nValid = 0;
  double l2NormVcm = 0.;
  double maxMagnitudeVcm = 0.;
};

struct ScreeningFieldChange {
  double relativeL2 = std::numeric_limits<double>::quiet_NaN();
  double relativeMax = std::numeric_limits<double>::quiet_NaN();
  std::size_t nCompared = 0;
};

ScreeningFieldGridSample SampleScreeningFieldGrid(
    Garfield::ComponentPoisson2d& sc, const SpaceChargeConfig& cfg);

ScreeningFieldChange CompareScreeningFieldGrids(
    const ScreeningFieldGridSample& previous,
    const ScreeningFieldGridSample& current);

ScreeningFieldSample SampleScreeningField(
    Garfield::ComponentPoisson2d& sc, const SpaceChargeConfig& cfg);

// Expect a few kV/cm; ~100 kV/cm means zExtentCm is far too small.
void ReportScreeningField(Garfield::ComponentPoisson2d& sc,
                          const SpaceChargeConfig& cfg);

// "lgad190V.sta" -> 190
double ParseBiasFromFilename(const std::string& path);

std::string FormatBias(double v);

bool SetImpactIonisationModel(Garfield::MediumSilicon& si,
                              const std::string& model);

void DumpTownsendCoefficients(
    Garfield::MediumSilicon& si, const std::string& csvPath,
    double eMinVcm = 1.e5, double eMaxVcm = 4.5e5,
    double eStepVcm = 2.5e4);

double MedianOf(std::vector<std::size_t> v);

// Stateless discard buffer; safe for concurrent OpenMP writes.
class NullStreamBuf : public std::streambuf {
 public:
  int overflow(int c) override { return c; }
  std::streamsize xsputn(const char*, std::streamsize n) override {
    return n;
  }
};

struct ReadoutStrip {
  std::string label;
  double centerUm, halfWidthUm;
};

bool StripsInsideMap(const std::vector<ReadoutStrip>& strips, double bx0,
                     double bx1);

void ScanValidity(Garfield::Component& cmp, double x0, double x1,
                  double y0, double y1, int nx, int ny,
                  const std::string& csvPath, double eps = 0.013e-4);

struct FieldProfile {
  double yTop = 0., yBot = 0., d = 0., yGain = 0., eMax = 0.;
  bool valid = false;
};

// 1D scan at x0: active-silicon bounds and gain-layer peak.
FieldProfile ScanFieldProfile(Garfield::Component& cmp, double x0,
                              double y0, double y1, int nScan,
                              const std::string& csvPath,
                              double eMinVcm = 100.);

struct FieldDumpResult {
  std::size_t nRows = 0, nBad = 0;
  double wMax = 0., wMaxX = 0., wMaxY = 0.;
};

// Columns: x_um,y_um,Ex_Vcm,Ey_Vcm,Emag_Vcm,V. Does not print.
FieldDumpResult DumpElectricField(Garfield::Component& cmp, double x0,
                                  double x1, double y0, double y1, int nx,
                                  int ny, const std::string& csvPath,
                                  double eps = 0.013e-4);

// Fine step in the gain-layer band, coarse elsewhere.
void ConfigureAvalanche(Garfield::AvalancheMC& av, double& stepCm,
                        double bulkStepCm, double yFineLoCm,
                        double yFineHiCm, double timeWindowNs,
                        std::size_t sizeCap, bool enableSignal,
                        bool enableMultithreading, const std::string& tag,
                        std::atomic<long long>& nCalls,
                        std::atomic<long long>& nFine,
                        std::atomic<long long>& nCoarse,
                        std::atomic<long long>& nPrinted);


void PrintWeightingSanity(Garfield::Component& wcmp,
                          const std::vector<ReadoutStrip>& strips, double yTop,
                          double yCollectHint, double yBot);

// Columns: x_um,y_um,<label>_phi,...
void DumpWeightingField(Garfield::Component& wcmp,
                        const std::vector<ReadoutStrip>& strips, double bx0,
                        double bx1, double yTop, double yBot,
                        const std::string& outPath, int nx = 250,
                        int ny = 150);

/* psi_k(x,y,t) = prompt + delayed, one CSV per time, plus
   <prefix>_summary.csv. Empty timesNs dumps the prompt map only. */
void DumpDynamicWeighting(Garfield::Component& wcmp,
                          const std::vector<ReadoutStrip>& strips, double bx0,
                          double bx1, double yLo, double yHi,
                          const std::vector<double>& timesNs,
                          const std::string& outPrefix, int nx = 200,
                          int ny = 120,
                          const std::vector<std::string>& extraLabels = {});

struct DelayedWeightingValidation {
  bool available = false;
  std::vector<double> timesNs;
  std::string message;
  double probedDelayedMax = 0.;  // 0 if probe skipped
};

/* Checks a shared, increasing delayed-weighting time grid AND a non-null
   delayed potential (DelayedSignalTimes alone can pass on an empty store).
   probeYCm > 0 enables the live probe. */
DelayedWeightingValidation ValidateDelayedWeightingData(
    Garfield::Component& weightingCmp,
    const std::vector<ReadoutStrip>& strips,
    double probeYCm = 0.);

// One Sensor + AvalancheMC per OpenMP thread; returns total electrons.
struct AvalanchePassConfig {
  // Optional drift-path CSV for animation, capped at driftMaxCarriers.
  std::ostream* driftCsv = nullptr;
  long driftMaxCarriers = 400;
  double fineStepCm = 5.e-6, bulkStepCm = 2.5e-5;
  double yFineLoCm = 0., yFineHiCm = 0.;
  double timeWindowNs = 4.;
  std::size_t sizeCap = 5000;
  bool multiplication = true;
  bool diffusion = true;
  bool enableSignal = true;
  double xMin = 0., yMin = 0., zMin = 0., xMax = 0., yMax = 0., zMax = 0.;
  double tStart = 0., tStep = 0.005;
  unsigned int nBins = 800;
  // Also hides genuine warnings like "not in a valid drift region".
  bool silenceGarfieldStderr = true;
  bool silenceGarfieldAvalancheStdout = true;
  bool silenceGarfieldSetupStdout = true;
  bool enableDelayedSignal = false;
  // Keep on for any AC-LGAD bipolar-signal run.
  bool requireDelayedWeightingData = true;
  std::vector<double> delayedSignalTimesNs;
  std::size_t delayedSignalAveragingOrder = 0;

  int ompDynamicChunk = 4;

  bool forceSerial = false;
};

unsigned long RunAvalanchePass(
    Garfield::Component& driftCmp,
    Garfield::Component& weightingCmp,
    // One per strip (ComponentGrid holds one electrode); empty = weightingCmp.
    const std::vector<Garfield::Component*>& perStripWeighting,
    const std::vector<ReadoutStrip>& strips,
    const std::vector<std::array<double, 4>>& primaries,
    const AvalanchePassConfig& pc,
    std::vector<std::vector<double>>& signalOut,
    std::ofstream* pairsCsv, unsigned long printEvery,
    const std::string& tag,
    // Space charge: scField adds to the field; scDeposit/scGrid collect charge.
    // Caller calls Solve().
    Garfield::Component* scField = nullptr,
    Garfield::ComponentPoisson2d* scDeposit = nullptr,
    const SpaceChargeConfig* scCfg = nullptr,
    ChargeGrid* scGrid = nullptr,
    std::vector<std::vector<double>>* promptSignalOut = nullptr,
    std::vector<std::vector<double>>* delayedSignalOut = nullptr,
    // Measures L_z; forces drift-line storage on.
    TransverseExtentAccumulator* lzAcc = nullptr,
    double lzBandLoCm = 0., double lzBandHiCm = 0.);

// G_e/G_eh vs step size; avalLadder must be configured via ConfigureAvalanche.
void RunConvergenceScan(Garfield::AvalancheMC& avalLadder, double& stepCm,
                        double x0, double yInj, std::size_t sizeCap,
                        const std::string& outDir);


struct FeedbackScanConfig {
  std::vector<std::string> models = {"okuto"};

  // Stops when SEM/mean < targetRelativeSem, within [minEvents, maxEvents].
  int minEvents = 500;
  int maxEvents = 5000;
  int batchSize = 250;
  double targetRelativeSem = 0.05;
  double heavyTailThresholdF = 3.0;

  bool runHoleSeed = true;
  bool runPairSeed = true;
};

/* Feedback scan modes: e_no_holes, e_full, h_full, eh_full.
   Caller must restore si's MIP model afterward. */
void RunModelComparison(
    Garfield::AvalancheMC& avalLadder, Garfield::MediumSilicon& si,
    double x0, double yInj, double stepCm, std::size_t ladderCap,
    const std::string& outDir, const std::string& biasLabel,
    double eMax, double gMax,
    const FeedbackScanConfig& cfg = FeedbackScanConfig{});
