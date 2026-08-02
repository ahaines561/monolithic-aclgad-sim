#pragma once

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <fstream>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "Garfield/AvalancheMC.hh"
#include "Garfield/Component.hh"
#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/ComponentPoisson2d.hh"
#include "Garfield/MediumSilicon.hh"

using Clock = std::chrono::steady_clock;

/* Silicon with a screening correction applied to impact ionisation only.

Garfield has no space-charge feedback: the avalanche never perturbs the
field it grows in. Silvaco's coupled Poisson solve does, and that is the
whole measured Garfield/Silvaco gain gap -- a ~2.4% (150 V) to ~3.1%
(180 V) effective deficit in the gain-layer field.

This scales the field seen by the Townsend calculation ONLY. Drift
velocity, diffusion and timing are untouched, which is deliberate: the
gain-OFF passes of the two codes already agree, so drift is validated and
must not be modified. AvalancheMC::GetTownsend falls through to the
medium when no Townsend map is loaded, and LoadSilvaco imports only the
potential and field, so this override is the multiplication path.

Set the scale as f = 1 - kappa * G, with kappa calibrated once against a
known Garfield/Silvaco pair (kappa = 0.00313 for lgad_deepjte, Okuto).
f = 1 restores stock behaviour exactly. */
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

/* AvalancheMC drifts carriers through a frozen field, so the avalanche never
screens the gain layer it grows in silvaco's coupled Poisson solve does,
and that is the measured ~2.4% (150 V) to ~3.1% (180 V) effective
gain-layer field deficit behind the Garfield/Silvaco gain gap.

ComponentPoisson2d solves Poisson on a triangular FEM mesh and accepts
deposited charge directly through AddCharge (elementary charges per cm of
depth) with a consistent finite-element load. It applies the medium's
relative permittivity properly, and Solve() reuses the existing
factorisation, so each feedback iteration costs one back-substitution
rather than a refactorisation.

Layering: Sensor::ElectricField sums over components, so

    sensor.AddComponent(&cmp);   // ComponentTcad2d, the Silvaco field
    sensor.AddComponent(&sc);    // ComponentPoisson2d, the perturbation

gives E_total = E_TCAD + E_screening. Add the TCAD component FIRST: the
sensor takes its medium from the first component returning status 0.

The TCAD map already carries the equilibrium space charge (ionised
dopants), so SetDoping is deliberately NOT called here. This component
must hold only the avalanche's excess carriers, or the depletion charge
is counted twice. */
struct SpaceChargeConfig {
  // FEM region outline [cm]; normally the field-map bounding box.
  // ComponentPoisson2d wants polygon edges parallel to x or y, so a
  // rectangle is the natural choice.
  double xMinCm = 0., xMaxCm = 0.;
  double yMinCm = 0., yMaxCm = 0.;
  double zMinCm = -5.e-4, zMaxCm = 5.e-4;

  // target element size near electrodes (hmin) and in the bulk (hmax) [cm]
  double hMinCm = 0.02e-4;
  double hMaxCm = 0.50e-4;

  /* Effective z-extent of the charge column [cm]. THIS IS THE PHYSICS
     INPUT, not a numerical knob. A 2D solve treats charge as uniform in
     z, and AddCharge takes elementary charges per cm of depth, so the
     carrier count is divided by this. Match it to the track
     cross-section Silvaco's SEU assumes; getting it wrong scales the
     whole screening effect linearly. */
  double zExtentCm = 1.e-4;

  // window over which the carrier density is time-averaged [ns]
  double tWindowNs = 1.;

  int maxIter = 6;           // max space-charge feedback iterations

  /* Under-relaxation:  rho_used = (1-lambda)*rho_prev + lambda*rho_new.
     The undamped map (lambda = 1) oscillates with contraction ratio
     r ~ -0.5 on this device, so it converges but slowly. The optimum is
     lambda = 1/(1-r) ~ 0.67; 0.5 is near-optimal with margin. */
  double relaxation = 0.5;

  /* Auxiliary grid for the spatial mixing, in nodes. The blend must be
     done on the charge DENSITY, not on a scalar amplitude: the cloud's
     shape changes between iterations (field/charge varied 44% across one
     measured run), so scaling the magnitude alone is not equivalent. */
  int mixNx = 101, mixNy = 501;

  /* The screening field is the fixed-point state variable. Convergence is
     checked on the full sampled screening field in the gain-layer
     neighbourhood. Gain remains an observable/diagnostic and is evaluated
     on the converged field in a final signal pass. Requiring three stable
     transitions avoids stopping on one accidental close pair. */
  double fieldTol = 0.01;
  int stableIterations = 3;
  int fieldSampleNx = 11;
  int fieldSampleNy = 41;
  double fieldSampleXHalfWidthCm = 5.e-4;  // 5 um
  double fieldSampleYHalfWidthCm = 1.e-4;  // 1 um

  /* Gain remains a diagnostic only; stochastic avalanche fluctuations do
     not gate field convergence. Relaxation provides noise suppression. */
  double gainTol = 0.01;     // diagnostic only; does not gate convergence

  // probe point for the screening-field diagnostic [cm]
  double xProbeCm = 0., yProbeCm = 0.;

  bool enabled = false;      // OFF by default: preserves stock behaviour
  bool verbose = true;
};

/* Build the FEM region, grounded electrodes and mesh, then Initialise().
   Doping is left unset (nd = na = 0) so the component carries only
   deposited charge. Returns false if Initialise() fails. */
bool SetupSpaceCharge(Garfield::ComponentPoisson2d& sc,
                      const SpaceChargeConfig& cfg,
                      Garfield::Medium* medium);

/* Deposit one avalanche's carriers, weighted by residence time.
   Electrons negative, holes positive; their separation is what opposes
   the gain-layer field. Does NOT call Solve() -- the caller does that
   once per iteration, after all avalanches have been deposited.

   REQUIRES aval.EnableDriftLines(true). AvalancheMC does not store drift
   paths otherwise (m_storeDriftLines defaults to false), so every
   EndPoint::path is empty and this silently deposits nothing.

   Returns the number of path segments that fell outside the mesh. */
std::size_t DepositAvalancheCharge(const Garfield::AvalancheMC& aval,
                                   Garfield::ComponentPoisson2d& sc,
                                   const SpaceChargeConfig& cfg);

/* Auxiliary grid holding the avalanche charge density so successive
   iterations can be blended spatially. Charge is deposited cloud-in-cell
   onto the four surrounding nodes, which conserves total charge and keeps
   sub-cell position information -- important because the gain layer is
   thin compared with any practical cell size. */
struct ChargeGrid {
  int nx = 0, ny = 0;
  double xMinCm = 0., xMaxCm = 0., yMinCm = 0., yMaxCm = 0.;
  // Cached inverse node spacing. Charge deposition is in the avalanche
  // hot path, so avoid recomputing two divisions for every drift segment.
  double invDx = 0., invDy = 0.;
  std::vector<double> rho;       // charge per node [e / cm of depth]

  void Init(const SpaceChargeConfig& cfg);
  void Clear() { std::fill(rho.begin(), rho.end(), 0.); }
  bool Contains(double xCm, double yCm) const noexcept {
    return xCm >= xMinCm && xCm <= xMaxCm &&
           yCm >= yMinCm && yCm <= yMaxCm;
  }
  void Deposit(double xCm, double yCm, double q);   // checked cloud-in-cell
  void DepositUnchecked(double xCm, double yCm, double q) noexcept;
  double Total() const;
  double AbsoluteTotal() const;
  double NodeX(int ix) const;
  double NodeY(int iy) const;
};

/* Accumulate one avalanche into the grid, residence-time weighted.
   Same physics as DepositAvalancheCharge but into the mixing grid
   instead of straight into the solver. Returns segments that fell
   outside the grid. */
std::size_t DepositAvalancheIntoGrid(const Garfield::AvalancheMC& aval,
                                     ChargeGrid& grid,
                                     const SpaceChargeConfig& cfg);

/* mixed = (1-lambda)*mixed + lambda*fresh, then push the result into the
   Poisson solver (ClearCharge + one AddCharge per non-empty node).
   Does NOT call Solve(). */
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

/* Sample the perturbation field on a fixed grid centred on the track and
   gain layer. This catches field-shape changes that a single probe cannot. */
ScreeningFieldGridSample SampleScreeningFieldGrid(
    Garfield::ComponentPoisson2d& sc, const SpaceChargeConfig& cfg);

/* Compare two samples made on the same grid. */
ScreeningFieldChange CompareScreeningFieldGrids(
    const ScreeningFieldGridSample& previous,
    const ScreeningFieldGridSample& current);

/* Sample the solved perturbation field and deposited line charge at the
   configured probe point. Useful for event-level CSV diagnostics. */
ScreeningFieldSample SampleScreeningField(
    Garfield::ComponentPoisson2d& sc, const SpaceChargeConfig& cfg);

/* Report the screening field at the probe point. Call after Solve().
   A few thousand V/cm against a ~376 kV/cm gain layer is the expected
   scale; ~100 kV/cm means zExtentCm is far too small. */
void ReportScreeningField(Garfield::ComponentPoisson2d& sc,
                          const SpaceChargeConfig& cfg);

// pulls e.g. 190 out of ".../lgad190V.sta" or "lgad_190V.sta"
double ParseBiasFromFilename(const std::string& path);

// "NA" if unparsed; whole numbers print without a decimal point
std::string FormatBias(double v);

/* Select one of the supported MediumSilicon impact-ionisation models.
   Returns false and prints an error for an unknown model name. */
bool SetImpactIonisationModel(Garfield::MediumSilicon& si,
                              const std::string& model);

/* Export the coefficients actually evaluated by Garfield at the current
   medium temperature/model. Field and coefficients are in V/cm and cm^-1. */
void DumpTownsendCoefficients(
    Garfield::MediumSilicon& si, const std::string& csvPath,
    double eMinVcm = 1.e5, double eMaxVcm = 4.5e5,
    double eStepVcm = 2.5e4);

double MedianOf(std::vector<std::size_t> v);

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
  bool valid = false;  // false if no valid drift medium found on the line
};

/* 1D scan at x0 from y0 to y1: finds the active-silicon bounds and the
gain-layer peak, writes the profile to a text file. Prints its own
diagnostics (active silicon bounds, peak field, low-field NOTE, or the
"no valid drift medium" error if !valid).*/

FieldProfile ScanFieldProfile(Garfield::Component& cmp, double x0,
                              double y0, double y1, int nScan,
                              const std::string& csvPath,
                              double eMinVcm = 100.);

struct FieldDumpResult {
  std::size_t nRows = 0, nBad = 0;
  double wMax = 0., wMaxX = 0., wMaxY = 0.;
};

/* E-field over an nx*ny grid spanning [x0,x1]x[y0,y1], to a text file:
x_um,y_um,Ex_Vcm,Ey_Vcm,Emag_Vcm,V. 

Does not print -- callers want different summary lines (window vs full-device); use the returned result to print. */
FieldDumpResult DumpElectricField(Garfield::Component& cmp, double x0,
                                  double x1, double y0, double y1, int nx,
                                  int ny, const std::string& csvPath,
                                  double eps = 0.013e-4);

/* Avalanche Stepping
- step function (fine band around gain layer, coarse elsewhere)
- size cap
- time window
- stepCm captured by ref so RunConvergenceScan can tune it */
void ConfigureAvalanche(Garfield::AvalancheMC& av, double& stepCm,
                        double bulkStepCm, double yFineLoCm,
                        double yFineHiCm, double timeWindowNs,
                        std::size_t sizeCap, bool enableSignal,
                        bool enableMultithreading, const std::string& tag,
                        std::atomic<long long>& nCalls,
                        std::atomic<long long>& nFine,
                        std::atomic<long long>& nCoarse,
                        std::atomic<long long>& nPrinted);


/* weighting field
- per-strip Ew/wpot startup check (probed just inside the gap) 
-plus an overlap sanity check under strip0. Pure diagnostic printing. */
void PrintWeightingSanity(Garfield::Component& wcmp,
                          const std::vector<ReadoutStrip>& strips, double yTop,
                          double yBot);

// per-strip weighting potential over a grid, to a text file:
// x_um,y_um,<label>_phi,...
void DumpWeightingField(Garfield::Component& wcmp,
                        const std::vector<ReadoutStrip>& strips, double bx0,
                        double bx1, double yTop, double yBot,
                        const std::string& outPath, int nx = 250,
                        int ny = 150);

struct DelayedWeightingValidation {
  bool available = false;
  std::vector<double> timesNs;
  std::string message;
  // Largest |delayed weighting potential| seen by the live probe. Stays 0 if
  // the probe was skipped (probeYCm <= 0).
  double probedDelayedMax = 0.;
};

/* Check that every requested electrode label has the same non-empty,
   strictly increasing delayed-weighting time grid, AND that the delayed
   weighting POTENTIAL is actually non-null there.

   The time-grid check alone is not sufficient: DelayedSignalTimes() falls back
   from the potential store to the field store, so it succeeds even when the
   store the signal path reads is empty. Pass probeYCm > 0 (a depth inside the
   bulk, in cm) to enable the live probe; pass 0 to skip it. */
DelayedWeightingValidation ValidateDelayedWeightingData(
    Garfield::Component& weightingCmp,
    const std::vector<ReadoutStrip>& strips,
    double probeYCm = 0.);

/* parallel avalanche pass (OpenMP)
- one Sensor + one AvalancheMC per thread; only the read-only field
  components are shared. Mirrors the pattern proven in 1000MIPs.cpp.
- per-thread signals are summed at the end (induced signals are additive)
- returns total electrons after multiplication
- compiles and runs correctly WITHOUT -fopenmp (pragmas ignored, serial) */
struct AvalanchePassConfig {
  double fineStepCm = 5.e-6, bulkStepCm = 2.5e-5;
  double yFineLoCm = 0., yFineHiCm = 0.;
  double timeWindowNs = 4.;
  std::size_t sizeCap = 5000;
  bool multiplication = true;
  // Diagnostic switch. AvalancheMC diffusion is enabled by default;
  // setting this false calls EnableDiffusion(false) for every local avalanche.
  bool diffusion = true;
  bool enableSignal = true;
  double xMin = 0., yMin = 0., zMin = 0., xMax = 0., yMax = 0., zMax = 0.;
  double tStart = 0., tStep = 0.005;
  unsigned int nBins = 800;
  // Garfield warns once per avalanche that it is running
  // single-threaded inside an OpenMP region (which is the intended
  // behaviour) -- hundreds of lines. Silencing also hides genuine
  // stderr warnings such as "not in a valid drift region".
  bool silenceGarfieldStderr = true;
  // Silence Garfield's per-avalanche transport progress meter and other
  // routine stdout emitted while the pass is running. Event-level summaries
  // are printed by main.cpp after RunAvalanchePass returns, so they remain
  // visible. Set false only for Garfield transport debugging.
  bool silenceGarfieldAvalancheStdout = true;
  // AddElectrode and SetTimeWindow write routine setup messages to stdout.
  // This narrower switch is used only when full-pass stdout silencing is off.
  bool silenceGarfieldSetupStdout = true;
  // Include the delayed weighting-field contribution in the signal.
  // The weighting component registered with Sensor::AddElectrode must provide
  // DelayedWeightingPotential/Field for each requested electrode label.
  bool enableDelayedSignal = false;
  // If true, an absent/inconsistent dynamic weighting map is an error.
  // Keep this enabled for any run advertised as an AC-LGAD bipolar-signal run.
  bool requireDelayedWeightingData = true;
  std::vector<double> delayedSignalTimesNs;
  std::size_t delayedSignalAveragingOrder = 0;

  // Dynamic OpenMP chunk. Small chunks improve load balance because
  // avalanche sizes vary strongly from one primary pair to another.
  int ompDynamicChunk = 4;

  // Optional serial execution for debugging or thread-safety studies.
  bool forceSerial = false;
};

unsigned long RunAvalanchePass(
    Garfield::Component& driftCmp,
    Garfield::Component& weightingCmp,
    // Optional: one weighting component PER STRIP, aligned with `strips`.
    // Required for ComponentGrid, which ignores the electrode label and can
    // therefore hold only one electrode per instance. When empty, every strip
    // is read from the single `weightingCmp`.
    const std::vector<Garfield::Component*>& perStripWeighting,
    const std::vector<ReadoutStrip>& strips,
    const std::vector<std::array<double, 4>>& primaries,
    const AvalanchePassConfig& pc,
    std::vector<std::vector<double>>& signalOut,
    std::ofstream* pairsCsv, unsigned long printEvery,
    const std::string& tag,
    /* Optional space-charge feedback. If scField is given it is added to
       every thread-local sensor, so carriers drift in E_TCAD + E_screen.
       If scDeposit and scCfg are given, each avalanche's carriers are
       deposited as they are produced (AddCharge is mutex-protected, so
       this is safe from inside the parallel loop). Solve() is NOT called
       here -- the caller does that between iterations. */
    Garfield::Component* scField = nullptr,
    Garfield::ComponentPoisson2d* scDeposit = nullptr,
    const SpaceChargeConfig* scCfg = nullptr,
    /* If scGrid is given, carriers are accumulated into the mixing grid
       instead of straight into the solver, so the caller can blend
       iterations. Each thread fills a private copy which is summed at
       the end, so this is safe in parallel. */
    ChargeGrid* scGrid = nullptr,
    /* Optional decomposition of the total signal. These arrays are filled
       only when non-null. delayedSignalOut is meaningful only when
       enableDelayedSignal=true and the weighting component supplies a
       dynamic weighting potential or field. */
    std::vector<std::vector<double>>* promptSignalOut = nullptr,
    std::vector<std::vector<double>>* delayedSignalOut = nullptr);

/* convergence scan
G_e/G_eh vs step size at a fixed injection point; avalLadder must
already be configured via ConfigureAvalanche. */
void RunConvergenceScan(Garfield::AvalancheMC& avalLadder, double& stepCm,
                        double x0, double yInj, std::size_t sizeCap,
                        const std::string& outDir);


struct FeedbackScanConfig {
  // Models to scan. Keep this to {"okuto"} for the primary
  // Garfield/Silvaco comparison; add vodm/massey/grant only as needed.
  std::vector<std::string> models = {"okuto"};

  // Adaptive stopping. Every mode runs at least minEvents and then stops
  // when SEM / mean reaches targetRelativeSem, or at maxEvents.
  int minEvents = 500;
  int maxEvents = 5000;
  int batchSize = 250;
  double targetRelativeSem = 0.05;
  double heavyTailThresholdF = 3.0;

  bool runHoleSeed = true;
  bool runPairSeed = true;
};

/* Impact-ionisation feedback scan.
   Modes:
   - e_no_holes: electron seed; generated holes counted but not transported
   - e_full: electron seed; generated holes transported and multiplied
   - h_full: hole seed; generated electrons transported
   - eh_full: primary electron-hole pair with full feedback

   Uses adaptive statistics, explicit heavy-tail/cap-limited statuses, and
   writes feedback_results_v2.csv, feedback_summary_v2.csv, and the raw
   per-mode avalanche-size text files.

   Caller must restore si's configured MIP model afterward if a later MIP
   pass is requested. */
void RunModelComparison(
    Garfield::AvalancheMC& avalLadder, Garfield::MediumSilicon& si,
    double x0, double yInj, double stepCm, std::size_t ladderCap,
    const std::string& outDir, const std::string& biasLabel,
    double eMax, double gMax,
    const FeedbackScanConfig& cfg = FeedbackScanConfig{});
