#pragma once

#include <atomic>
#include <chrono>
#include <string>
#include <vector>

#include "Garfield/AvalancheMC.hh"
#include "Garfield/Component.hh"
#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/MediumSilicon.hh"

using Clock = std::chrono::steady_clock;


double ElapsedS(const Clock::time_point& t0);

// pulls e.g. 190 out of ".../lgad190V.sta" or "lgad_190V.sta"
double ParseBiasFromFilename(const std::string& path);

// "NA" if unparsed; whole numbers print without a decimal point
std::string FormatBias(double v);

double MedianOf(std::vector<std::size_t> v);

struct Strip {
  std::string label;
  double centerUm, halfWidthUm;
};

bool StripsInsideMap(const std::vector<Strip>& strips, double bx0,
                     double bx1);

void ScanValidity(Garfield::Component& cmp, double x0, double x1,
                  double y0, double y1, int nx, int ny,
                  const std::string& csvPath);

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
                              const std::string& csvPath);

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
void PrintWeightingSanity(Garfield::ComponentAnalyticField& wcmp,
                          const std::vector<Strip>& strips, double yTop,
                          double yBot);

// per-strip weighting potential over a grid, to a text file:
// x_um,y_um,<label>_phi,...
void DumpWeightingField(Garfield::ComponentAnalyticField& wcmp,
                        const std::vector<Strip>& strips, double bx0,
                        double bx1, double yTop, double yBot,
                        const std::string& outPath, int nx = 250,
                        int ny = 150);

/* convergence scan
G_e/G_eh vs step size at a fixed injection point; avalLadder must
already be configured via ConfigureAvalanche. */
void RunConvergenceScan(Garfield::AvalancheMC& avalLadder, double& stepCm,
                        double x0, double yInj, std::size_t sizeCap,
                        const std::string& outDir);


/* impact-ionisation model comparison
- G_e/G_eh for vodm/okuto/massey/grant at one step size and injection point
- writes eh/ge_sizes*.txt and appends to results.csv
- Caller must reset si's model afterward. */
void RunModelComparison(Garfield::AvalancheMC& avalLadder,
                        Garfield::MediumSilicon& si, double x0,
                        double yInj, double stepCm, std::size_t ladderCap,
                        const std::string& outDir,
                        const std::string& biasLabel, double eMax,
                        double gMax);
