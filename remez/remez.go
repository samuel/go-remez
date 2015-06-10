package remez

import (
	"errors"
	"math"
)

/********************************************************
 *
 *  Code taken from remez.c by Erik Kvaleberg which was
 *    converted from an original FORTRAN by
 *
 * AUTHORS: JAMES H. MCCLELLAN
 *
 *         DEPARTMENT OF ELECTRICAL ENGINEERING AND COMPUTER SCIENCE
 *         MASSACHUSETTS INSTITUTE OF TECHNOLOGY
 *         CAMBRIDGE, MASS. 02139
 *
 *         THOMAS W. PARKS
 *         DEPARTMENT OF ELECTRICAL ENGINEERING
 *         RICE UNIVERSITY
 *         HOUSTON, TEXAS 77001
 *
 *         LAWRENCE R. RABINER
 *         BELL LABORATORIES
 *         MURRAY HILL, NEW JERSEY 07974
 *
 *
 *  Adaptation to C by
 *      egil kvaleberg
 *      husebybakken 14a
 *      0379 oslo, norway
 *  Email:
 *      egil@kvaleberg.no
 *  Web:
 *      http://www.kvaleberg.com/
 *
 *
 *********************************************************/

type FilterType int

const (
	BandPass FilterType = iota
	Differentiator
	Hilbert
)

const (
	pi  = math.Pi
	pi2 = math.Pi * 2
)

/*
 *-----------------------------------------------------------------------
 * FUNCTION: lagrange_interp (d)
 *  FUNCTION TO CALCULATE THE LAGRANGE INTERPOLATION
 *  COEFFICIENTS FOR USE IN THE FUNCTION gee.
 *-----------------------------------------------------------------------
 */
func lagrangeInterp(k, n, m int, x []float64) float64 {
	retval := 1.0
	q := x[k]

	for l := 1; l <= m; l++ {
		for j := l; j <= n; j += m {
			if j != k {
				retval *= 2.0 * (q - x[j])
			}
		}
	}
	return 1.0 / retval
}

/*
 *-----------------------------------------------------------------------
 * FUNCTION: freq_eval (gee)
 *  FUNCTION TO EVALUATE THE FREQUENCY RESPONSE USING THE
 *  LAGRANGE INTERPOLATION FORMULA IN THE BARYCENTRIC FORM
 *-----------------------------------------------------------------------
 */
func freqEval(k, n int, grid, x, y, ad []float64) float64 {
	d := 0.0
	p := 0.0
	xf := math.Cos(pi2 * grid[k])

	for j := 1; j <= n; j++ {
		c := ad[j] / (xf - x[j])
		d += c
		p += c * y[j]
	}

	return p / d
}

/*
 *-----------------------------------------------------------------------
 * SUBROUTINE: remez
 *  THIS SUBROUTINE IMPLEMENTS THE REMEZ EXCHANGE ALGORITHM
 *  FOR THE WEIGHTED CHEBYSHEV APPROXIMATION OF A CONTINUOUS
 *  FUNCTION WITH A SUM OF COSINES.  INPUTS TO THE SUBROUTINE
 *  ARE A DENSE GRID WHICH REPLACES THE FREQUENCY AXIS, THE
 *  DESIRED FUNCTION ON THIS GRID, THE WEIGHT FUNCTION ON THE
 *  GRID, THE NUMBER OF COSINES, AND AN INITIAL GUESS OF THE
 *  EXTREMAL FREQUENCIES.  THE PROGRAM MINIMIZES THE CHEBYSHEV
 *  ERROR BY DETERMINING THE BSMINEST LOCATION OF THE EXTREMAL
 *  FREQUENCIES (POINTS OF MAXIMUM ERROR) AND THEN CALCULATES
 *  THE COEFFICIENTS OF THE BEST APPROXIMATION.
 *-----------------------------------------------------------------------
 */
func remez(des, grid, bands, wt []float64, ngrid int, iext []int, alpha []float64, nfcns, itrmax int, dimsize int) (float64, error) {
	a := make([]float64, dimsize+1)
	p := make([]float64, dimsize+1)
	q := make([]float64, dimsize+1)
	ad := make([]float64, dimsize+1)
	x := make([]float64, dimsize+1)
	y := make([]float64, dimsize+1)

	devl := -1.0
	nz := nfcns + 1
	nzz := nfcns + 2

	var comp, dev, y1 float64

IterationLoop:
	for niter := 0; niter <= itrmax; niter++ {
		if niter == itrmax {
			return dev, errors.New("remez: reached max iterations")
		}

		iext[nzz] = ngrid + 1

		for j := 1; j <= nz; j++ {
			x[j] = math.Cos(grid[iext[j]] * pi2)
		}

		jet := (nfcns-1)/15 + 1
		for j := 1; j <= nz; j++ {
			ad[j] = lagrangeInterp(j, nz, jet, x)
		}

		dnum, dden := 0.0, 0.0
		for j, k := 1, 1.0; j <= nz; j, k = j+1, -k {
			l := iext[j]
			dnum += ad[j] * des[l]
			dden += k * ad[j] / wt[l]
		}
		dev = dnum / dden

		/* printf("DEVIATION = %lg\n",*dev); */

		nu := 1.0
		if dev > 0.0 {
			nu = -1.0
		}
		dev = math.Abs(dev) // dev = -nu * dev
		for j, k := 1, nu; j <= nz; j, k = j+1, -k {
			l := iext[j]
			y[j] = des[l] + k*dev/wt[l]
		}
		if dev <= devl {
			/* finished */
			return dev, errors.New("remez: deviation decreased")
		}
		devl = dev

		jchnge := 0
		k1 := iext[1]
		knz := iext[nz]
		klow := 0
		nut := -nu

		down := func(l, j int) {
			for {
				l--
				if l <= klow {
					break
				}
				e := (freqEval(l, nz, grid, x, y, ad) - des[l]) * wt[l]
				if nut*e-comp <= 0.0 {
					break
				}
				comp = nut * e
			}
			klow = iext[j]
			iext[j] = l + 1
			jchnge++
		}

		up := func(l, j, kup int) {
			for {
				l++
				if l >= kup {
					break
				}
				e := (freqEval(l, nz, grid, x, y, ad) - des[l]) * wt[l]
				if nut*e-comp <= 0.0 {
					break
				}
				comp = nut * e
			}
			iext[j] = l - 1
			klow = l - 1
			jchnge++
		}

		/*
		 * SEARCH FOR THE EXTREMAL FREQUENCIES OF THE BEST APPROXIMATION
		 */

		for j := 1; j < nzz; j++ {
			kup := iext[j+1]
			nut = -nut
			if j == 2 {
				y1 = comp
			}
			comp = dev

			l := iext[j] + 1
			if l < kup {
				e := (freqEval(l, nz, grid, x, y, ad) - des[l]) * wt[l]
				if nut*e-comp > 0.0 {
					comp = nut * e
					up(l, j, kup)
					continue
				}
			}

			l--

			for {
				l--
				if l <= klow {
					l = iext[j] + 1
					if jchnge > 0 {
						iext[j] = l - 1
						klow = l - 1
						jchnge++
					} else {
						for {
							l++
							if l >= kup {
								klow = iext[j]
								break
							}
							e := (freqEval(l, nz, grid, x, y, ad) - des[l]) * wt[l]
							if nut*e-comp > 0.0 {
								comp = nut * e
								up(l, j, kup)
								break
							}
						}
					}
					break
				}
				e := (freqEval(l, nz, grid, x, y, ad) - des[l]) * wt[l]
				if nut*e-comp > 0.0 {
					comp = nut * e
					down(l, j)
					break
				}
				if jchnge > 0 {
					klow = iext[j]
					break
				}
			}
		}

		if k1 > iext[1] {
			k1 = iext[1]
		}
		if knz < iext[nz] {
			knz = iext[nz]
		}

		luck := 6
		nut1 := nut
		nut = -nu
		comp *= 1.00001
		j := nzz

		for l := 1; l < k1; l++ {
			e := (freqEval(l, nz, grid, x, y, ad) - des[l]) * wt[l]
			if nut*e-comp > 0.0 {
				comp = nut * e
				up(l, j, k1)
				j = nzz + 1
				luck = 1

				if comp > y1 {
					y1 = comp
				}
				k1 = iext[nzz]
				break
			}
		}

		klow = knz
		nut = -nut1
		comp = y1 * 1.00001

		for l := ngrid; l > klow; l-- {
			e := (freqEval(l, nz, grid, x, y, ad) - des[l]) * wt[l]
			if nut*e-comp > 0.0 {
				comp = nut * e
				down(l, j)

				kn := iext[nzz]
				for i := 1; i <= nfcns; i++ {
					iext[i] = iext[i+1]
				}
				iext[nz] = kn
				continue IterationLoop
			}
		}

		if luck != 6 {
			for i := 1; i <= nfcns; i++ {
				iext[nzz-i] = iext[nz-i]
			}
			iext[1] = k1
		} else if jchnge <= 0 {
			break
		}
	}

	/*
	 *    CALCULATION OF THE COEFFICIENTS OF THE BEST APPROXIMATION
	 *    USING THE INVERSE DISCRETE FOURIER TRANSFORM
	 */
	nm1 := nfcns - 1
	fsh := 1.0e-06
	gtemp := grid[1]
	x[nzz] = -2.0
	cn := float64(2*nfcns - 1)
	delf := 1.0 / cn
	l := 1
	kkk := 0

	if bands[0] == 0.0 && bands[len(bands)-1] == 0.5 {
		kkk = 1
	}

	if nfcns <= 3 {
		kkk = 1
	}

	var aa, bb float64
	if kkk != 1 {
		dtemp := math.Cos(pi2 * grid[1])
		dnum := math.Cos(pi2 * grid[ngrid])
		aa = 2.0 / (dtemp - dnum)
		bb = -(dtemp + dnum) / (dtemp - dnum)
	}

	for j := 1; j <= nfcns; j++ {
		ft := float64(j-1) * delf
		xt := math.Cos(pi2 * ft)
		if kkk != 1 {
			xt = (xt - bb) / aa
			// /*XX* ckeck up !! */
			// xt1 = sqrt(1.0-xt*xt);
			// ft = atan2(xt1,xt)/pi2;

			ft = math.Acos(xt) / pi2
		}
		for {
			xe := x[l]
			if xt > xe {
				if (xt - xe) < fsh {
					a[j] = y[l]
					break
				}
				grid[1] = ft
				a[j] = freqEval(1, nz, grid, x, y, ad)
				break
			}
			if (xe - xt) < fsh {
				a[j] = y[l]
				break
			}
			l++
		}
		if l > 1 {
			l = l - 1
		}
	}

	grid[1] = gtemp
	dden := pi2 / cn
	for j := 1; j <= nfcns; j++ {
		dtemp := 0.0
		dnum := float64(j-1) * dden
		if nm1 >= 1 {
			for k := 1; k <= nm1; k++ {
				dtemp += a[k+1] * math.Cos(dnum*float64(k))
			}
		}
		alpha[j] = 2.0*dtemp + a[1]
	}

	for j := 2; j <= nfcns; j++ {
		alpha[j] *= 2.0 / cn
	}
	alpha[1] /= cn

	if kkk != 1 {
		p[1] = 2.0*alpha[nfcns]*bb + alpha[nm1]
		p[2] = 2.0 * aa * alpha[nfcns]
		q[1] = alpha[nfcns-2] - alpha[nfcns]
		for j := 2; j <= nm1; j++ {
			if j >= nm1 {
				aa *= 0.5
				bb *= 0.5
			}
			p[j+1] = 0.0
			for k := 1; k <= j; k++ {
				a[k] = p[k]
				p[k] = 2.0 * bb * a[k]
			}
			p[2] += a[1] * 2.0 * aa
			for k := 1; k <= j-1; k++ {
				p[k] += q[k] + aa*a[k+1]
			}
			for k := 3; k <= j+1; k++ {
				p[k] += aa * a[k-1]
			}

			if j != nm1 {
				for k := 1; k <= j; k++ {
					q[k] = -a[k]
				}
				q[1] += alpha[nfcns-1-j]
			}
		}
		for j := 1; j <= nfcns; j++ {
			alpha[j] = p[j]
		}
	}

	if nfcns <= 3 {
		alpha[nfcns+1] = 0.0
		alpha[nfcns+2] = 0.0
	}
	return dev, nil
}

/*
 *-----------------------------------------------------------------------
 * FUNCTION: eff
 *  FUNCTION TO CALCULATE THE DESIRED MAGNITUDE RESPONSE
 *  AS A FUNCTION OF FREQUENCY.
 *  AN ARBITRARY FUNCTION OF FREQUENCY CAN BE
 *  APPROXIMATED IF THE USER REPLACES THIS FUNCTION
 *  WITH THE APPROPRIATE CODE TO EVALUATE THE IDEAL
 *  MAGNITUDE.  NOTE THAT THE PARAMETER FREQ IS THE
 *  VALUE OF NORMALIZED FREQUENCY NEEDED FOR EVALUATION.
 *-----------------------------------------------------------------------
 */
func eff(freq float64, fx []float64, lband int, filterType FilterType) float64 {
	if filterType != Differentiator {
		return fx[lband]
	}
	return fx[lband] * freq
}

/*
 *-----------------------------------------------------------------------
 * FUNCTION: wate
 *  FUNCTION TO CALCULATE THE WEIGHT FUNCTION AS A FUNCTION
 *  OF FREQUENCY.  SIMILAR TO THE FUNCTION eff, THIS FUNCTION CAN
 *  BE REPLACED BY A USER-WRITTEN ROUTINE TO CALCULATE ANY
 *  DESIRED WEIGHTING FUNCTION.
 *-----------------------------------------------------------------------
 */
func wate(freq float64, fx, wtx []float64, lband int, filterType FilterType) float64 {
	if filterType != Differentiator {
		return wtx[lband]
	}
	if fx[lband] >= 0.0001 {
		return wtx[lband] / freq
	}
	return wtx[lband]
}

func Remez(numtaps int, bands, response, weight []float64, filterType FilterType, maxiter, gridDensity int) ([]float64, error) {
	var change float64
	var ngrid int

	dimsize := int(math.Ceil(float64(numtaps)/2.0 + 2))
	wrksize := gridDensity * dimsize
	nbands := len(bands) / 2
	/* Note:  code assumes these arrays start at 1 */
	h := make([]float64, numtaps+1)
	des := make([]float64, wrksize+1)
	grid := make([]float64, wrksize+1)
	wt := make([]float64, wrksize+1)
	alpha := make([]float64, dimsize+1)
	iext := make([]int, dimsize+1)

	/* Set up problem on dense grid */

	neg := true
	if filterType == BandPass {
		neg = false
	}
	nodd := numtaps%2 == 1
	nfcns := numtaps / 2
	if nodd && !neg {
		nfcns++
	}

	/*
	 * SET UP THE DENSE GRID. THE NUMBER OF POINTS IN THE GRID
	 * IS (FILTER LENGTH + 1)*GRID DENSITY/2
	 */
	grid[1] = bands[0]
	delf := 0.5 / float64(gridDensity*nfcns)
	if neg {
		if bands[0] < delf {
			grid[1] = delf
		}
	}

	/*
	 * CALCULATE THE DESIRED MAGNITUDE RESPONSE AND THE WEIGHT
	 * FUNCTION ON THE GRID
	 */
	j := 1
	l := 1
	lband := 0
	for {
		fup := bands[l]
		for {
			temp := grid[j]
			des[j] = eff(temp, response, lband, filterType)
			wt[j] = wate(temp, response, weight, lband, filterType)
			j++
			if j > wrksize {
				return nil, errors.New("remez: too many points or too dense grid")
			}
			grid[j] = temp + delf
			if grid[j] > fup {
				break
			}
		}

		grid[j-1] = fup
		des[j-1] = eff(fup, response, lband, filterType)
		wt[j-1] = wate(fup, response, weight, lband, filterType)
		lband++
		l += 2
		if lband >= nbands {
			break
		}
		grid[j] = bands[l-1]
	}

	ngrid = j - 1
	if neg == nodd {
		if grid[ngrid] > 0.5-delf {
			ngrid--
		}
	}

	grid = grid[:ngrid+1]

	/*
	 * SET UP A NEW APPROXIMATION PROBLEM WHICH IS EQUIVALENT
	 * TO THE ORIGINAL PROBLEM
	 */
	if !neg {
		if !nodd {
			for j := 1; j <= ngrid; j++ {
				change = math.Cos(pi * grid[j])
				des[j] /= change
				wt[j] *= change
			}
		}
	} else {
		if nodd {
			for j := 1; j <= ngrid; j++ {
				change = math.Sin(pi2 * grid[j])
				des[j] /= change
				wt[j] *= change
			}
		} else {
			for j := 1; j <= ngrid; j++ {
				change = math.Sin(pi * grid[j])
				des[j] /= change
				wt[j] *= change
			}
		}
	}

	/*XX*/
	temp := float64(ngrid-1) / float64(nfcns)
	for j := 1; j <= nfcns; j++ {
		iext[j] = int((float64(j-1) * temp) + 1) /* round? !! */
	}
	iext[nfcns+1] = ngrid
	nm1 := nfcns - 1
	nz := nfcns + 1

	if _, err := remez(des, grid, bands, wt, ngrid, iext, alpha, nfcns, maxiter, dimsize); err != nil {
		return nil, err
	}

	/*
	 * CALCULATE THE IMPULSE RESPONSE.
	 */
	if !neg {
		if nodd {
			for j := 1; j <= nm1; j++ {
				h[j] = 0.5 * alpha[nz-j]
			}
			h[nfcns] = alpha[1]
		} else {
			h[1] = 0.25 * alpha[nfcns]
			for j := 2; j <= nm1; j++ {
				h[j] = 0.25 * (alpha[nz-j] + alpha[nfcns+2-j])
			}
			h[nfcns] = 0.5*alpha[1] + 0.25*alpha[2]
		}
	} else {
		if nodd {
			h[1] = 0.25 * alpha[nfcns]
			h[2] = 0.25 * alpha[nm1]
			for j := 3; j <= nm1; j++ {
				h[j] = 0.25 * (alpha[nz-j] - alpha[nfcns+3-j])
			}
			h[nfcns] = 0.5*alpha[1] - 0.25*alpha[3]
			h[nz] = 0.0
		} else {
			h[1] = 0.25 * alpha[nfcns]
			for j := 2; j <= nm1; j++ {
				h[j] = 0.25 * (alpha[nz-j] - alpha[nfcns+2-j])
			}
			h[nfcns] = 0.5*alpha[1] - 0.25*alpha[2]
		}
	}

	for j := 1; j <= nfcns; j++ {
		k := numtaps + 1 - j
		if !neg {
			h[k] = h[j]
		} else {
			h[k] = -h[j]
		}
	}
	if neg && nodd {
		h[nz] = 0.0
	}

	return h[1:], nil
}
