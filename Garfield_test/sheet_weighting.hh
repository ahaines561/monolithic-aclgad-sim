#ifndef SHEET_WEIGHTING_HH
#define SHEET_WEIGHTING_HH

// Dynamic weighting potentials for a strip AC-LGAD (header-only).
// Unit voltage step on one electrode at t = 0:
//
//   bulk :  Laplace(phi) = 0                          0 <= y <= d
//   sheet:  d(sigma)/dt = (1/R_sq) d2u/dx2            y = 0
//           sigma = eps_Si * (DtN u) + C_ox(x) * (u - V_pad(x))
//
// The slab's DtN map is diagonal in cos(k_m x), eigenvalue k_m coth(k_m d),
// so phi_m(y) = u_m sinh(k_m (d - y)) / sinh(k_m d) and only the sheet is
// solved. t -> 0+: sigma = 0 (n+ transparent); t -> inf: psi -> 0.
// Crank-Nicolson on a log time grid.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace mass {

// dense partial-pivot LU
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

// Lengths in cm; defaults from aclgad.sta.
struct SheetConfig {
  double xMin = 0.;
  double xMax = 5.0e-2;
  double ySheet = 0.430466e-4;   // depletion edge
  double yBack = 5.0e-3;
  double tOx = 1.0e-5;

  double epsSi = 11.9;
  double epsOx = 3.9;            // SiO2
  double rSheet = 1500.;         // ohm/sq

  // same order as the readout strips
  std::vector<std::pair<double, double>> pads = {
      {30.0e-4, 70.0e-4}, {220.0e-4, 270.0e-4}, {430.0e-4, 470.0e-4}};
  // ohmic DC contact
  double dcX0 = 0.;
  double dcX1 = 3.0e-4;

  // Grounded metal on oxide (DC electrode, capacitively coupled).
  std::vector<std::pair<double, double>> groundedPads = {{3.0e-4, 5.0e-4}};

  std::size_t nSheet = 301;      // sheet nodes = cosine modes
  std::size_t nOutY = 81;
  std::size_t stepsPerDecade = 60;
  double tStartNs = 1.0e-4;

  std::vector<double> timesNs = {
      0.002, 0.005, 0.01, 0.02, 0.035, 0.05, 0.075, 0.1, 0.2, 0.35, 0.5,
      0.75, 1., 1.5, 2., 3., 5., 7.5, 10., 15., 20., 30., 50., 75.,
      100., 150.};
};

class SheetWeighting {
 public:
  explicit SheetWeighting(const SheetConfig& cfg) : m_cfg(cfg) { Build(); }

  std::size_t NumElectrodes() const { return m_cfg.pads.size(); }
  const std::vector<double>& TimesNs() const { return m_cfg.timesNs; }

  // psi_k(x, y, 0+)
  double Prompt(std::size_t k, double x, double y) const {
    return Sample(m_psi[k], 0, x, y);
  }

  // psi_k(t) - psi_k(0+): Garfield adds the prompt term itself.
  double Delayed(std::size_t k, double x, double y, double tNs) const {
    const auto& ts = m_cfg.timesNs;
    if (tNs <= ts.front()) {
      const double f = tNs <= 0. ? 0. : tNs / ts.front();
      return f * (Sample(m_psi[k], 1, x, y) - Sample(m_psi[k], 0, x, y));
    }
    if (tNs >= ts.back())
      return Sample(m_psi[k], ts.size(), x, y) - Sample(m_psi[k], 0, x, y);
    // binary search keeps cost independent of the slice count
    std::size_t i = static_cast<std::size_t>(
        std::lower_bound(ts.begin(), ts.end(), tNs) - ts.begin());
    i = std::max<std::size_t>(1, std::min(i, ts.size() - 1));
    const double f = (tNs - ts[i - 1]) / (ts[i] - ts[i - 1]);
    const double a = Sample(m_psi[k], i, x, y);
    const double b = Sample(m_psi[k], i + 1, x, y);
    return (1. - f) * a + f * b - Sample(m_psi[k], 0, x, y);
  }

  /* Completeness (sum psi = 1) and decay checks. Both pass for ANY R_sq or
     t_ox, so they cannot validate sheet parameters; use
     tools/compare_silvaco_waveform.py for that. */
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

    // C_ox weighted by cell coverage so pad widths are grid-independent.
    m_cOxVal = cOxVal;
    auto coverage = [&](std::size_t i, double a, double b) {
      const double lo = std::max(m_cfg.xMin, m_x[i] - 0.5 * m_dx);
      const double hi = std::min(m_cfg.xMax, m_x[i] + 0.5 * m_dx);
      const double w = hi - lo;
      if (w <= 0.) return 0.;
      return std::max(0., std::min(hi, b) - std::max(lo, a)) / w;
    };

    m_cOx.assign(n, 0.);
    m_padFrac.assign(m_cfg.pads.size(), std::vector<double>(n, 0.));
    for (std::size_t k = 0; k < m_cfg.pads.size(); ++k)
      for (std::size_t i = 0; i < n; ++i) {
        const double f = coverage(i, m_cfg.pads[k].first, m_cfg.pads[k].second);
        m_padFrac[k][i] = f;
        m_cOx[i] += cOxVal * f;
      }

    // grounded metal: loads every pad, driven only in the DC solve
    m_dcFrac.assign(n, 0.);
    for (const auto& g : m_cfg.groundedPads)
      for (std::size_t i = 0; i < n; ++i) {
        const double f = coverage(i, g.first, g.second);
        m_dcFrac[i] += f;
        m_cOx[i] += cOxVal * f;
      }

    // Dirichlet has no fractional form: pin cells > half covered.
    m_isDC.assign(n, 0);
    for (std::size_t i = 0; i < n; ++i)
      if (coverage(i, m_cfg.dcX0, m_cfg.dcX1) > 0.5) m_isDC[i] = 1;

    m_k.resize(n);
    m_decay.assign(n * m_cfg.nOutY, 0.);
    std::vector<double> C(n * n);
    for (std::size_t m = 0; m < n; ++m) {
      m_k[m] = M_PI * static_cast<double>(m) / L;
      for (std::size_t i = 0; i < n; ++i)
        C[m * n + i] = std::cos(m_k[m] * (m_x[i] - m_cfg.xMin));
      for (std::size_t j = 0; j < m_cfg.nOutY; ++j) {
        const double yy = m_y[j] - m_cfg.ySheet;
        m_decay[m * m_cfg.nOutY + j] = ModeProfile(m_k[m], d, yy);
      }
    }
    // Cinv: nodal -> modal
    std::vector<double> Ct(n * n);
    for (std::size_t i = 0; i < n; ++i)
      for (std::size_t m = 0; m < n; ++m) Ct[i * n + m] = C[m * n + i];
    m_Cinv = Invert(Ct, n);
    m_C = std::move(Ct);   // now nodal[i] = sum_m C[i*n+m] * modal[m]

    // DtN eigenvalues k coth(kd) -> 1/d as k -> 0
    std::vector<double> lam(n);
    for (std::size_t m = 0; m < n; ++m) {
      const double kd = m_k[m] * d;
      lam[m] = (kd < 1e-8) ? 1.0 / d
                           : m_k[m] * std::cosh(std::min(kd, 700.)) /
                                 std::sinh(std::min(kd, 700.));
      if (kd > 700.) lam[m] = m_k[m];
    }
    // B = epsSi * C * diag(lam) * Cinv
    m_B.assign(n * n, 0.);
    for (std::size_t i = 0; i < n; ++i)
      for (std::size_t m = 0; m < n; ++m) {
        const double w = epsSi * m_C[i * n + m] * lam[m];
        if (w == 0.) continue;
        for (std::size_t j = 0; j < n; ++j)
          m_B[i * n + j] += w * m_Cinv[m * n + j];
      }

    // sheet conductance, reflecting edges
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
    // electrodes: pads..., DC contact, backplane
    const std::size_t nEl = m_cfg.pads.size() + 2;
    m_psi.assign(nEl, std::vector<double>((m_cfg.timesNs.size() + 1) *
                                          m_cfg.nOutY * n, 0.));

    std::vector<double> grid;
    BuildTimeGrid(grid);

    for (std::size_t k = 0; k < nEl; ++k) {
      const bool isDCel = (k == m_cfg.pads.size());
      const bool isBP = (k == m_cfg.pads.size() + 1);
      // drive = C_ox of the cell fraction under this electrode
      std::vector<double> cOxDrive(n, 0.);
      if (isDCel) {
        for (std::size_t i = 0; i < n; ++i)
          cOxDrive[i] = m_cOxVal * m_dcFrac[i];
      } else if (!isBP) {
        for (std::size_t i = 0; i < n; ++i)
          cOxDrive[i] = m_cOxVal * m_padFrac[k][i];
      }
      const double vDC = isDCel ? 1. : 0.;
      const double vBP = isBP ? 1. : 0.;

      // t = 0+: sigma = 0  ->  (B + C_ox) u = C_ox V_pad + eps V_bp / d
      std::vector<double> M(m_B);
      std::vector<double> rhs(n);
      const double dSlab = m_cfg.yBack - m_cfg.ySheet;
      const double eps0b = 8.8541878128e-14;
      for (std::size_t i = 0; i < n; ++i) {
        M[i * n + i] += m_cOx[i];
        rhs[i] = cOxDrive[i] + eps0b * m_cfg.epsSi * vBP / dSlab;
      }
      ApplyDirichlet(M, rhs, vDC, n);
      std::vector<double> u = LU(M, n).Solve(rhs);
      Store(k, 0, u, vBP);

      // Crank-Nicolson on d(sigma)/dt = A u (constant drive cancels)
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

  // psi(x, y) from the sheet potential + backplane linear term
  void Store(std::size_t k, std::size_t it, const std::vector<double>& u,
             double vBP) {
    const std::size_t n = m_cfg.nSheet;
    const double d = m_cfg.yBack - m_cfg.ySheet;
    const std::vector<double>& uu = u;
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
  std::vector<std::vector<double>> m_padFrac;  // [pad][node] coverage
  std::vector<double> m_dcFrac;
  double m_cOxVal = 0.;
  std::vector<char> m_isDC;
  std::vector<double> m_B, m_A, m_C, m_Cinv, m_decay;
  std::vector<std::vector<double>> m_psi;
};

}  // namespace mass

#endif  // SHEET_WEIGHTING_HH
