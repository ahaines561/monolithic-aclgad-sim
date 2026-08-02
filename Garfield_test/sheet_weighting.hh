#ifndef SHEET_WEIGHTING_HH
#define SHEET_WEIGHTING_HH

// Dynamic weighting potentials for a strip AC-LGAD, computed in-process.
//
// Header-only and dependency-free: no Eigen, no external solver, no data
// files, no build-system changes. Include it and construct one object.
// For a unit voltage step on one electrode at t = 0 we solve
//
//   bulk :  Laplace(phi) = 0        on 0 <= y <= d   (depleted, uniform eps)
//   sheet:  d(sigma)/dt = (1/R_sq) d2u/dx2           at y = 0
//
// where u(x, t) is the potential of the undepleted n+ and sigma its free
// surface charge density,
//
//   sigma = eps_Si * (DtN u)  +  C_ox(x) * (u - V_pad(x)).
//
// The bulk is a homogeneous slab with reflecting boundaries at x = 0 and
// x = L and a grounded backplane at y = d, so its Dirichlet-to-Neumann map is
// DIAGONAL in the cosine basis cos(k_m x), k_m = m*pi/L, with eigenvalues
// k_m * coth(k_m d)  ->  1/d as k -> 0. That removes the bulk mesh entirely:
// the y dependence is analytic,
//
//   phi_m(y) = u_m * sinh(k_m (d - y)) / sinh(k_m d),
//
// and the whole problem collapses to a dense n_x system on the sheet.
//
// Two limits fix the physics and are checked by SelfTest():
//   t -> 0+   no lateral charge has moved, so sigma = 0 on every free sheet
//             node. That is a Neumann condition: the n+ is electrostatically
//             TRANSPARENT and does NOT sit at ground in the gaps.
//   t -> inf  the sheet is tied through R to the grounded DC contact, so
//             u -> 0, psi -> 0 everywhere, and the induced charge on every
//             AC pad integrates to zero. That is the bipolar signal, and it
//             follows from the boundary conditions rather than being imposed.
//
// Time integration is Crank-Nicolson on an internal logarithmic grid.
// Backward Euler damps stiff modes only algebraically and grossly
// overestimates the tail when the step is large; CN is second order and was
// verified to converge from the opposite side to the same limit.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace mass {

// minimal dense linear algebra (partial-pivot LU)
class LU {
 public:
  LU() = default;
  explicit LU(std::vector<double> a, std::size_t n) : m_n(n), m_a(std::move(a)) {
    m_piv.resize(n);
    for (std::size_t i = 0; i < n; ++i) m_piv[i] = i;
    for (std::size_t k = 0; k < n; ++k) {
      std::size_t p = k;
      double big = std::abs(m_a[k * n + k]);
      for (std::size_t i = k + 1; i < n; ++i) {
        const double v = std::abs(m_a[i * n + k]);
        if (v > big) { big = v; p = i; }
      }
      if (big < 1e-300) throw std::runtime_error("LU: singular matrix");
      if (p != k) {
        for (std::size_t j = 0; j < n; ++j)
          std::swap(m_a[k * n + j], m_a[p * n + j]);
        std::swap(m_piv[k], m_piv[p]);
      }
      const double d = m_a[k * n + k];
      for (std::size_t i = k + 1; i < n; ++i) {
        const double f = m_a[i * n + k] / d;
        m_a[i * n + k] = f;
        if (f == 0.) continue;
        for (std::size_t j = k + 1; j < n; ++j)
          m_a[i * n + j] -= f * m_a[k * n + j];
      }
    }
  }

  std::vector<double> Solve(const std::vector<double>& b) const {
    const std::size_t n = m_n;
    std::vector<double> x(n);
    for (std::size_t i = 0; i < n; ++i) x[i] = b[m_piv[i]];
    for (std::size_t i = 1; i < n; ++i)
      for (std::size_t j = 0; j < i; ++j) x[i] -= m_a[i * n + j] * x[j];
    for (std::size_t i = n; i-- > 0;) {
      for (std::size_t j = i + 1; j < n; ++j) x[i] -= m_a[i * n + j] * x[j];
      x[i] /= m_a[i * n + i];
    }
    return x;
  }

 private:
  std::size_t m_n = 0;
  std::vector<double> m_a;
  std::vector<std::size_t> m_piv;
};

// configuration
struct SheetConfig {
  // Geometry, all in cm. Defaults are the values extracted from aclgad.sta.
  double xMin = 0.;
  double xMax = 5.0e-2;          // 500 um
  double ySheet = 0.430466e-4;   // depletion edge: top of the depleted bulk
  double yBack = 5.0e-3;         // backplane, 50 um
  double tOx = 1.0e-5;           // coupling dielectric, 100 nm

  double epsSi = 11.9;
  double epsOx = 3.9;            // 3.9 = SiO2, 7.5 = Si3N4. VERIFY THIS.
  double rSheet = 1500.;         // ohm / square

  // Pad footprints in cm, in the same order as the readout strips.
  std::vector<std::pair<double, double>> pads = {
      {30.0e-4, 70.0e-4}, {220.0e-4, 270.0e-4}, {430.0e-4, 470.0e-4}};
  // Grounded DC contact on the n+ (the `bulk` electrode in the .sta).
  double dcX0 = 0.;
  double dcX1 = 5.0e-4;

  std::size_t nSheet = 301;      // sheet nodes (and cosine modes)
  std::size_t nOutY = 81;        // output grid in y
  std::size_t stepsPerDecade = 60;
  double tStartNs = 1.0e-4;

  std::vector<double> timesNs = {
      0.002, 0.005, 0.01, 0.02, 0.035, 0.05, 0.075, 0.1, 0.2, 0.35, 0.5,
      0.75, 1., 1.5, 2., 3., 5., 7.5, 10., 15., 20., 30., 50., 75.,
      100., 150.};
};

// solver
class SheetWeighting {
 public:
  explicit SheetWeighting(const SheetConfig& cfg) : m_cfg(cfg) { Build(); }

  std::size_t NumElectrodes() const { return m_cfg.pads.size(); }
  const std::vector<double>& TimesNs() const { return m_cfg.timesNs; }

  // Prompt weighting potential, psi_k(x, y, 0+).
  double Prompt(std::size_t k, double x, double y) const {
    return Sample(m_psi[k], 0, x, y);
  }

  // Delayed REMAINDER, psi_k(x, y, t) - psi_k(x, y, 0+). This is the quantity
  // Garfield's delayed machinery expects: Sensor adds the prompt term and the
  // delayed term, so storing psi(t) itself would double-count the prompt and
  // break I_total = I_prompt + I_delayed.
  double Delayed(std::size_t k, double x, double y, double tNs) const {
    const auto& ts = m_cfg.timesNs;
    if (tNs <= ts.front()) {
      const double f = tNs <= 0. ? 0. : tNs / ts.front();
      return f * (Sample(m_psi[k], 1, x, y) - Sample(m_psi[k], 0, x, y));
    }
    if (tNs >= ts.back())
      return Sample(m_psi[k], ts.size(), x, y) - Sample(m_psi[k], 0, x, y);
    std::size_t i = 1;
    while (i + 1 < ts.size() && ts[i] < tNs) ++i;
    const double f = (tNs - ts[i - 1]) / (ts[i] - ts[i - 1]);
    const double a = Sample(m_psi[k], i, x, y);
    const double b = Sample(m_psi[k], i + 1, x, y);
    return (1. - f) * a + f * b - Sample(m_psi[k], 0, x, y);
  }

  // Completeness: stepping every electrode together drives no current, so
  // phi == 1 and the weighting potentials of ALL electrodes (pads, DC
  // contact, backplane) must sum to exactly 1 at every (x, y, t). This one
  // check validates the DtN map, the sheet balance, the Dirichlet handling
  // and the time stepping simultaneously.
  bool SelfTest(std::ostream& os = std::cout) const {
    bool ok = true;
    const std::size_t nt = m_cfg.timesNs.size() + 1;
    double worst = 0.;
    for (std::size_t it = 0; it < nt; ++it) {
      for (std::size_t j = 0; j < m_cfg.nOutY; j += 8) {
        for (std::size_t i = 0; i < m_cfg.nSheet; i += 8) {
          double s = 0.;
          for (std::size_t k = 0; k < m_psi.size(); ++k)
            s += m_psi[k][it * m_cfg.nOutY * m_cfg.nSheet +
                          j * m_cfg.nSheet + i];
          worst = std::max(worst, std::abs(s - 1.));
        }
      }
    }
    os << "  [sheet] completeness  max|sum psi - 1| = " << worst
       << (worst < 1e-6 ? "   PASS" : "   FAIL") << "\n";
    ok &= worst < 1e-6;

    for (std::size_t k = 0; k < m_cfg.pads.size(); ++k) {
      const double p0 = MaxAbs(m_psi[k], 0);
      const double pT = MaxAbs(m_psi[k], m_cfg.timesNs.size());
      os << "  [sheet] electrode " << k << "  max|psi(0+)| = " << p0
         << "   max|psi(t_end)|/max|psi(0+)| = " << pT / p0
         << (pT / p0 < 0.02 ? "   PASS" : "   FAIL") << "\n";
      ok &= pT / p0 < 0.02;
    }
    return ok;
  }

 private:
  //construction 
  void Build() {
    const std::size_t n = m_cfg.nSheet;
    const double L = m_cfg.xMax - m_cfg.xMin;
    const double d = m_cfg.yBack - m_cfg.ySheet;
    m_dx = L / static_cast<double>(n - 1);
    m_x.resize(n);
    for (std::size_t i = 0; i < n; ++i)
      m_x[i] = m_cfg.xMin + static_cast<double>(i) * m_dx;
    m_y.resize(m_cfg.nOutY);
    for (std::size_t j = 0; j < m_cfg.nOutY; ++j)
      m_y[j] = m_cfg.ySheet + d * static_cast<double>(j) /
                                  static_cast<double>(m_cfg.nOutY - 1);

    const double eps0 = 8.8541878128e-14;   // F/cm
    const double epsSi = eps0 * m_cfg.epsSi;
    const double cOxVal = eps0 * m_cfg.epsOx / m_cfg.tOx;

    // Per-node oxide capacitance and pad membership.
    m_cOx.assign(n, 0.);
    m_padOf.assign(n, -1);
    for (std::size_t k = 0; k < m_cfg.pads.size(); ++k)
      for (std::size_t i = 0; i < n; ++i)
        if (m_x[i] >= m_cfg.pads[k].first && m_x[i] <= m_cfg.pads[k].second) {
          m_cOx[i] = cOxVal;
          m_padOf[i] = static_cast<int>(k);
        }
    m_isDC.assign(n, 0);
    for (std::size_t i = 0; i < n; ++i)
      if (m_x[i] >= m_cfg.dcX0 && m_x[i] <= m_cfg.dcX1) m_isDC[i] = 1;

    // Cosine basis: reflecting in x at both ends, so cos(k_m x) with
    // k_m = m*pi/L diagonalises the bulk.
    m_k.resize(n);
    m_decay.assign(n * m_cfg.nOutY, 0.);
    std::vector<double> C(n * n);
    for (std::size_t m = 0; m < n; ++m) {
      m_k[m] = M_PI * static_cast<double>(m) / L;
      for (std::size_t i = 0; i < n; ++i)
        C[m * n + i] = std::cos(m_k[m] * (m_x[i] - m_cfg.xMin));
      // Analytic y profile of this mode, sinh(k(d-y))/sinh(kd).
      for (std::size_t j = 0; j < m_cfg.nOutY; ++j) {
        const double yy = m_y[j] - m_cfg.ySheet;
        m_decay[m * m_cfg.nOutY + j] = ModeProfile(m_k[m], d, yy);
      }
    }
    // Cinv maps nodal values -> modal amplitudes.
    std::vector<double> Ct(n * n);
    for (std::size_t i = 0; i < n; ++i)
      for (std::size_t m = 0; m < n; ++m) Ct[i * n + m] = C[m * n + i];
    m_Cinv = Invert(Ct, n);
    m_C = std::move(Ct);   // now nodal[i] = sum_m C[i*n+m] * modal[m]

    // Dirichlet-to-Neumann: diagonal eigenvalue k coth(kd), -> 1/d as k -> 0.
    std::vector<double> lam(n);
    for (std::size_t m = 0; m < n; ++m) {
      const double kd = m_k[m] * d;
      lam[m] = (kd < 1e-8) ? 1.0 / d
                           : m_k[m] * std::cosh(std::min(kd, 700.)) /
                                 std::sinh(std::min(kd, 700.));
      if (kd > 700.) lam[m] = m_k[m];
    }
    // B = epsSi * C * diag(lam) * Cinv, dense.
    m_B.assign(n * n, 0.);
    for (std::size_t i = 0; i < n; ++i)
      for (std::size_t m = 0; m < n; ++m) {
        const double w = epsSi * m_C[i * n + m] * lam[m];
        if (w == 0.) continue;
        for (std::size_t j = 0; j < n; ++j)
          m_B[i * n + j] += w * m_Cinv[m * n + j];
      }

    // Sheet conductance operator, conservative form, reflecting at the edges.
    m_A.assign(n * n, 0.);
    const double g = 1.0 / (m_cfg.rSheet * m_dx * m_dx);
    for (std::size_t i = 0; i < n; ++i) {
      const std::size_t im = (i > 0) ? i - 1 : 1;
      const std::size_t ip = (i + 1 < n) ? i + 1 : n - 2;
      m_A[i * n + im] += g;
      m_A[i * n + ip] += g;
      m_A[i * n + i] -= 2. * g;
    }

    Solve();
  }

  static double ModeProfile(double k, double d, double y) {
    if (k * d < 1e-8) return 1.0 - y / d;
    if (k * d > 700.) return std::exp(-k * y);   // avoid sinh overflow
    return std::sinh(k * (d - y)) / std::sinh(k * d);
  }

  static std::vector<double> Invert(const std::vector<double>& a,
                                    std::size_t n) {
    LU lu(a, n);
    std::vector<double> inv(n * n);
    std::vector<double> e(n, 0.);
    for (std::size_t c = 0; c < n; ++c) {
      std::fill(e.begin(), e.end(), 0.);
      e[c] = 1.;
      const auto col = lu.Solve(e);
      for (std::size_t r = 0; r < n; ++r) inv[r * n + c] = col[r];
    }
    return inv;
  }

  void Solve() {
    const std::size_t n = m_cfg.nSheet;
    // Electrodes: pads, then the DC contact, then the backplane. The last two
    // are solved only so SelfTest can check completeness.
    const std::size_t nEl = m_cfg.pads.size() + 2;
    m_psi.assign(nEl, std::vector<double>((m_cfg.timesNs.size() + 1) *
                                          m_cfg.nOutY * n, 0.));

    std::vector<double> grid;
    BuildTimeGrid(grid);

    for (std::size_t k = 0; k < nEl; ++k) {
      const bool isDCel = (k == m_cfg.pads.size());
      const bool isBP = (k == m_cfg.pads.size() + 1);
      std::vector<double> vPad(n, 0.);
      if (!isDCel && !isBP)
        for (std::size_t i = 0; i < n; ++i)
          if (m_padOf[i] == static_cast<int>(k)) vPad[i] = 1.;
      const double vDC = isDCel ? 1. : 0.;
      const double vBP = isBP ? 1. : 0.;

      // t = 0+ : sigma = 0 on every free node  ->  (B + C_ox) u = C_ox V_pad,
      // with the backplane offset folded in as a uniform shift.
      std::vector<double> M(m_B);
      std::vector<double> rhs(n);
      const double dSlab = m_cfg.yBack - m_cfg.ySheet;
      const double eps0b = 8.8541878128e-14;
      for (std::size_t i = 0; i < n; ++i) {
        M[i * n + i] += m_cOx[i];
        // sigma = eps*(DtN u) - eps*V_bp/d + C_ox*(u - V_pad); setting
        // sigma = 0 puts both drives on the right-hand side.
        rhs[i] = m_cOx[i] * vPad[i] + eps0b * m_cfg.epsSi * vBP / dSlab;
      }
      ApplyDirichlet(M, rhs, vDC, n);
      std::vector<double> u = LU(M, n).Solve(rhs);
      Store(k, 0, u, vBP);

      // Crank-Nicolson: d(sigma)/dt = A u, sigma = (B + C_ox) u - C_ox V_pad,
      // and the drive is constant for t > 0 so it cancels in the difference.
      std::size_t out = 1;
      double tPrev = 0.;
      for (double t : grid) {
        const double dt = (t - tPrev) * 1e-9;
        tPrev = t;
        {
          std::vector<double> lhs(n * n), r(n, 0.);
          for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j) {
              const double s = m_B[i * n + j] + (i == j ? m_cOx[i] : 0.);
              lhs[i * n + j] = s / dt - 0.5 * m_A[i * n + j];
              r[i] += (s / dt + 0.5 * m_A[i * n + j]) * u[j];
            }
          ApplyDirichlet(lhs, r, vDC, n);
          u = LU(lhs, n).Solve(r);
        }
        if (out <= m_cfg.timesNs.size() &&
            std::abs(t - m_cfg.timesNs[out - 1]) < 1e-12 * std::max(1., t)) {
          Store(k, out, u, vBP);
          ++out;
        }
      }
      if (out != m_cfg.timesNs.size() + 1)
        throw std::runtime_error("SheetWeighting: internal time grid missed "
                                 "a requested output time");
    }
  }

  void BuildTimeGrid(std::vector<double>& grid) const {
    const double tmax = m_cfg.timesNs.back();
    const double nd = std::log10(tmax / m_cfg.tStartNs);
    const std::size_t np = static_cast<std::size_t>(
        std::ceil(nd * static_cast<double>(m_cfg.stepsPerDecade)));
    for (std::size_t i = 0; i <= np; ++i)
      grid.push_back(m_cfg.tStartNs *
                     std::pow(10., nd * static_cast<double>(i) /
                                       static_cast<double>(np)));
    for (double t : m_cfg.timesNs) grid.push_back(t);
    std::sort(grid.begin(), grid.end());
    grid.erase(std::unique(grid.begin(), grid.end(),
                           [](double a, double b) {
                             return std::abs(a - b) < 1e-12 * std::max(1., a);
                           }),
               grid.end());
  }

  void ApplyDirichlet(std::vector<double>& M, std::vector<double>& rhs,
                      double vDC, std::size_t n) const {
    for (std::size_t i = 0; i < n; ++i) {
      if (!m_isDC[i]) continue;
      for (std::size_t j = 0; j < n; ++j) M[i * n + j] = 0.;
      M[i * n + i] = 1.;
      rhs[i] = vDC;
    }
  }

  // Reconstruct psi(x, y) from the sheet potential using the analytic mode
  // profiles, and add the backplane's linear contribution.
  void Store(std::size_t k, std::size_t it, const std::vector<double>& u,
             double vBP) {
    const std::size_t n = m_cfg.nSheet;
    const double d = m_cfg.yBack - m_cfg.ySheet;
    const std::vector<double>& uu = u;   // w(sheet) = u by construction
    std::vector<double> mod(n, 0.);
    for (std::size_t m = 0; m < n; ++m) {
      double s = 0.;
      for (std::size_t j = 0; j < n; ++j) s += m_Cinv[m * n + j] * uu[j];
      mod[m] = s;
    }
    auto* dst = &m_psi[k][it * m_cfg.nOutY * n];
    for (std::size_t j = 0; j < m_cfg.nOutY; ++j) {
      const double yy = m_y[j] - m_cfg.ySheet;
      for (std::size_t i = 0; i < n; ++i) {
        double s = 0.;
        for (std::size_t m = 0; m < n; ++m)
          s += mod[m] * m_C[i * n + m] * m_decay[m * m_cfg.nOutY + j];
        dst[j * n + i] = s + vBP * (yy / d);
      }
    }
  }

  double MaxAbs(const std::vector<double>& v, std::size_t it) const {
    const std::size_t sz = m_cfg.nOutY * m_cfg.nSheet;
    double r = 0.;
    for (std::size_t i = 0; i < sz; ++i)
      r = std::max(r, std::abs(v[it * sz + i]));
    return r;
  }

  double Sample(const std::vector<double>& v, std::size_t it, double x,
                double y) const {
    const std::size_t n = m_cfg.nSheet;
    const double fx = (x - m_cfg.xMin) / m_dx;
    const double dy = (m_cfg.yBack - m_cfg.ySheet) /
                      static_cast<double>(m_cfg.nOutY - 1);
    const double fy = (y - m_cfg.ySheet) / dy;
    if (fx < 0. || fy < 0. || fx > static_cast<double>(n - 1) ||
        fy > static_cast<double>(m_cfg.nOutY - 1))
      return 0.;
    const auto i = static_cast<std::size_t>(fx);
    const auto j = static_cast<std::size_t>(fy);
    const std::size_t i1 = std::min(i + 1, n - 1);
    const std::size_t j1 = std::min(j + 1, m_cfg.nOutY - 1);
    const double a = fx - static_cast<double>(i);
    const double b = fy - static_cast<double>(j);
    const auto* p = &v[it * m_cfg.nOutY * n];
    return (1 - a) * (1 - b) * p[j * n + i] + a * (1 - b) * p[j * n + i1] +
           (1 - a) * b * p[j1 * n + i] + a * b * p[j1 * n + i1];
  }

  SheetConfig m_cfg;
  double m_dx = 0.;
  std::vector<double> m_x, m_y, m_k, m_cOx;
  std::vector<int> m_padOf;
  std::vector<char> m_isDC;
  std::vector<double> m_B, m_A, m_C, m_Cinv, m_decay;
  std::vector<std::vector<double>> m_psi;
};

}  // namespace mass

#endif  // SHEET_WEIGHTING_HH
