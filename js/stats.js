/* EpiTables — statistics core (no dependencies) */

const Stats = (() => {

  /* --- distribution functions --- */

  function normalCDF(x) {
    const t = 1 / (1 + 0.2316419 * Math.abs(x));
    const d = 0.3989423 * Math.exp(-x * x / 2);
    const p = d * t * (0.3193815 + t * (-0.3565638 + t * (1.7814779 + t * (-1.8212560 + t * 1.3302744))));
    return x > 0 ? 1 - p : p;
  }

  function chiSqCDF(x, k) {
    if (x <= 0) return 0;
    return regularizedGammaP(k / 2, x / 2);
  }

  function chiSqP(chi, df) {
    return 1 - chiSqCDF(chi, df);
  }

  function betaCDF(x, a, b) {
    if (x <= 0) return 0;
    if (x >= 1) return 1;
    return regularizedIncompleteBeta(x, a, b);
  }

  function betaInv(p, a, b) {
    // Newton-Raphson on regularized incomplete beta
    let x = a / (a + b);
    for (let i = 0; i < 200; i++) {
      const fx = betaCDF(x, a, b) - p;
      const dfx = Math.exp((a - 1) * Math.log(x) + (b - 1) * Math.log(1 - x) - logBeta(a, b));
      const dx = fx / dfx;
      x -= dx;
      x = Math.max(1e-10, Math.min(1 - 1e-10, x));
      if (Math.abs(dx) < 1e-10) break;
    }
    return x;
  }

  function regularizedGammaP(a, x) {
    if (x < 0) return 0;
    if (x === 0) return 0;
    if (x < a + 1) {
      // series
      let sum = 1 / a, term = 1 / a;
      for (let n = 1; n < 300; n++) {
        term *= x / (a + n);
        sum += term;
        if (Math.abs(term) < Math.abs(sum) * 1e-12) break;
      }
      return sum * Math.exp(-x + a * Math.log(x) - logGamma(a));
    } else {
      // continued fraction
      let b = x + 1 - a, c = 1e30, d = 1 / b, h = d;
      for (let i = 1; i <= 300; i++) {
        const an = -i * (i - a);
        b += 2;
        d = an * d + b; if (Math.abs(d) < 1e-30) d = 1e-30;
        c = b + an / c; if (Math.abs(c) < 1e-30) c = 1e-30;
        d = 1 / d;
        const del = d * c;
        h *= del;
        if (Math.abs(del - 1) < 1e-12) break;
      }
      return 1 - Math.exp(-x + a * Math.log(x) - logGamma(a)) * h;
    }
  }

  function regularizedIncompleteBeta(x, a, b) {
    if (x < (a + 1) / (a + b + 2)) {
      return betaCF(x, a, b) * Math.exp(a * Math.log(x) + b * Math.log(1 - x) - logBeta(a, b)) / a;
    } else {
      return 1 - betaCF(1 - x, b, a) * Math.exp(b * Math.log(1 - x) + a * Math.log(x) - logBeta(a, b)) / b;
    }
  }

  function betaCF(x, a, b) {
    const qab = a + b, qap = a + 1, qam = a - 1;
    let c = 1, d = 1 - qab * x / qap;
    if (Math.abs(d) < 1e-30) d = 1e-30;
    d = 1 / d;
    let h = d;
    for (let m = 1; m <= 200; m++) {
      const m2 = 2 * m;
      let aa = m * (b - m) * x / ((qam + m2) * (a + m2));
      d = 1 + aa * d; if (Math.abs(d) < 1e-30) d = 1e-30;
      c = 1 + aa / c; if (Math.abs(c) < 1e-30) c = 1e-30;
      d = 1 / d; h *= d * c;
      aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
      d = 1 + aa * d; if (Math.abs(d) < 1e-30) d = 1e-30;
      c = 1 + aa / c; if (Math.abs(c) < 1e-30) c = 1e-30;
      d = 1 / d;
      const del = d * c; h *= del;
      if (Math.abs(del - 1) < 1e-12) break;
    }
    return h;
  }

  function logGamma(x) {
    const c = [76.18009172947146,-86.50532032941677,24.01409824083091,-1.231739572450155,0.1208650973866179e-2,-0.5395239384953e-5];
    let y = x, tmp = x + 5.5;
    tmp -= (x + 0.5) * Math.log(tmp);
    let ser = 1.000000000190015;
    for (let j = 0; j < 6; j++) { y++; ser += c[j] / y; }
    return -tmp + Math.log(2.5066282746310005 * ser / x);
  }

  function logBeta(a, b) { return logGamma(a) + logGamma(b) - logGamma(a + b); }

  function logFact(n) {
    let s = 0;
    for (let i = 2; i <= n; i++) s += Math.log(i);
    return s;
  }

  /* --- CI methods --- */

  function wilson95(x, n) {
    if (n === 0) return [0, 0];
    const z = 1.96, p = x / n;
    const denom = 1 + z * z / n;
    const center = (p + z * z / (2 * n)) / denom;
    const margin = z * Math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom;
    return [Math.max(0, center - margin), Math.min(1, center + margin)];
  }

  function wald95(x, n) {
    if (n === 0) return [0, 1];
    const p = x / n, se = Math.sqrt(p * (1 - p) / n);
    return [Math.max(0, p - 1.96 * se), Math.min(1, p + 1.96 * se)];
  }

  function clopperpearson95(x, n) {
    const lo = x === 0 ? 0 : betaInv(0.025, x, n - x + 1);
    const hi = x === n ? 1 : betaInv(0.975, x + 1, n - x);
    return [lo, hi];
  }

  function logCI95(est, a, b, c, d) {
    if (a === 0 || b === 0 || c === 0 || d === 0) return [null, null];
    const se = Math.sqrt(1/a + 1/b + 1/c + 1/d);
    return [est * Math.exp(-1.96 * se), est * Math.exp(1.96 * se)];
  }

  function rrLogCI95(a, b, c, d) {
    const n1 = a + b, n2 = c + d;
    if (n1 === 0 || n2 === 0 || a === 0 || c === 0) return [null, null];
    const rr = (a / n1) / (c / n2);
    const se = Math.sqrt(b / (a * n1) + d / (c * n2));
    return [rr * Math.exp(-1.96 * se), rr * Math.exp(1.96 * se)];
  }

  function rdCI95(a, b, c, d) {
    const n1 = a + b, n2 = c + d;
    if (n1 === 0 || n2 === 0) return [null, null];
    const p1 = a / n1, p2 = c / n2, rd = p1 - p2;
    const se = Math.sqrt(p1 * (1 - p1) / n1 + p2 * (1 - p2) / n2);
    return [rd - 1.96 * se, rd + 1.96 * se];
  }

  /* --- tests --- */

  function chiSquare(a, b, c, d) {
    const N = a + b + c + d;
    if (N === 0) return { chi2: null, p: null };
    const chi2 = N * Math.pow(a * d - b * c, 2) / ((a+b) * (c+d) * (a+c) * (b+d));
    return { chi2, p: chiSqP(chi2, 1) };
  }

  function fisherExact(a, b, c, d) {
    const logP = (a, b, c, d) =>
      logFact(a+b) + logFact(c+d) + logFact(a+c) + logFact(b+d)
      - logFact(a) - logFact(b) - logFact(c) - logFact(d) - logFact(a+b+c+d);
    const obs = Math.exp(logP(a, b, c, d));
    let p = 0;
    const r1 = a+b, col1 = a+c, n = a+b+c+d;
    for (let aa = 0; aa <= Math.min(r1, col1); aa++) {
      const bb = r1-aa, cc = col1-aa, dd = n-r1-cc;
      if (bb < 0 || cc < 0 || dd < 0) continue;
      const prob = Math.exp(logP(aa, bb, cc, dd));
      if (prob <= obs + 1e-10) p += prob;
    }
    return Math.min(p, 1);
  }

  function mcnemar(b, c, n) {
    if (b + c === 0) return null;
    const chi2_corr = Math.pow(Math.abs(b - c) - 1, 2) / (b + c);
    const chi2_uncorr = Math.pow(b - c, 2) / (b + c);
    const se_diff = Math.sqrt((b + c - Math.pow(b - c, 2) / n)) / n;
    return {
      chi2_corr, p_corr: chiSqP(chi2_corr, 1),
      chi2_uncorr, p_uncorr: chiSqP(chi2_uncorr, 1),
      se_diff
    };
  }

  /* --- formatting --- */

  function fmt(x, dec = 2) {
    if (x === null || x === undefined || isNaN(x) || !isFinite(x)) return '—';
    return x.toFixed(dec);
  }

  function fmtP(p) {
    if (p === null || isNaN(p)) return '—';
    if (p < 0.001) return '<0.001';
    return p.toFixed(3);
  }

  function pct(x, dec = 1) {
    if (x === null || isNaN(x) || !isFinite(x)) return '—';
    return (x * 100).toFixed(dec) + '%';
  }

  function fmtCI(lo, hi, dec = 2) {
    if (lo === null || hi === null) return '';
    return `95% CI: ${fmt(lo, dec)} – ${fmt(hi, dec)}`;
  }

  return {
    normalCDF, chiSqP, betaInv,
    wilson95, wald95, clopperpearson95,
    logCI95, rrLogCI95, rdCI95,
    chiSquare, fisherExact, mcnemar,
    fmt, fmtP, pct, fmtCI
  };
})();
