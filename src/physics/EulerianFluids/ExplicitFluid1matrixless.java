/**
 * Port of Benedikt Bitterli incremental fluids source available at,
 * https://github.com/tunabrain/incremental-fluids
 *
 * Copyright (c) 2013 Benedikt Bitterli
 *
 * This software is provided 'as-is', without any express or implied warranty.
 * In no event will the authors be held liable for any damages arising from the
 * use of this software.
 *
 * Permission is granted to anyone to use this software for any purpose,
 * including commercial applications, and to alter it and redistribute it
 * freely, subject to the following restrictions:
 *
 * 1. The origin of this software must not be misrepresented; you must not claim
 * that you wrote the original software. If you use this software in a product,
 * an acknowledgment in the product documentation would be appreciated but is
 * not required.
 *
 * 2. Altered source versions must be plainly marked as such, and must not be
 * misrepresented as being the original software.
 *
 * 3. This notice may not be removed or altered from any source distribution.
 */
package physics.EulerianFluids;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.sqrt;
import static java.lang.String.format;
import static java.util.logging.Level.ALL;
import static java.util.logging.Level.INFO;
import static java.util.logging.Level.SEVERE;
import java.util.logging.Logger;
import static java.util.logging.Logger.getLogger;
import static math.Interpolation.lerp;
import org.apache.commons.math3.geometry.euclidean.twod.Vector2D;

/**
 * @author Gonçalo Amador {@literal &} Abel Gomes
 */
class ExplicitFluid1matrixless {

    static final Logger log = getLogger("ExplicitFluid1matrixless");

    static {
        log.setLevel(ALL);
    }

    public static void main(String args[]) {
        /**
         * Play with these constants, if you want
         */

        int sizeX = 128;
        int sizeY = 128;

        double density = 0.1;
        double timestep = 0.005;

        FluidSolver solver = new FluidSolver(sizeX, sizeY, density);

        double time = 0.0;
        int iterations = 0;

        while (time < 8.0) {
            /**
             * Use four substeps per iteration
             */
            for (int i = 0; i < 4; i++) {
                solver.addInflow(0.45, 0.2, 0.1, 0.01, 1.0, 0.0, 3.0);
                solver.update(timestep);
                time += timestep;
            }

            iterations++;
            File f = new File("Frame" + iterations + ".ppm");
            try {
                solver.toImage(f);
            } catch (IOException ex) {
                log.log(SEVERE, null, ex);
            }
        }
    }

    /**
     * This is the class representing fluid quantities such as density and
     * velocity on the MAC grid. It saves attributes such as offset from the top
     * left grid cell, grid width and height as well as cell size.
     *
     * It also contains two memory buffers: A source (_src) buffer and a
     * destination (_dst) buffer. Most operations on fluid quantities can be
     * done in-place; that is, they write to the same buffer they're reading
     * from (which is always _src). However, some operations, such as advection,
     * cannot be done in-place. Instead, they will write to the _dst buffer.
     * Once the operation is completed, flip() can be called to swap the source
     * and destination buffers, such that the result of the operation is visible
     * to subsequent operations.
     */
    static class FluidQuantity {

        /**
         * Memory buffers for fluid quantity
         */
        double[] _src;
        double[] _dst;
        double[] _temp;

        /**
         * Width and height
         */
        int _w;
        int _h;

        /**
         * X and Y offset from top left grid cell. This is (0.5,0.5) for
         * centered quantities such as density, and (0.0, 0.5) or (0.5, 0.0) for
         * jittered quantities like the velocity.
         */
        double _ox;
        double _oy;
        /**
         * Grid cell size
         */
        double _hx;

        FluidQuantity(int w, int h, double ox, double oy, double hx) {
            _w = w;
            _h = h;
            _ox = ox;
            _oy = oy;
            _hx = hx;

            _src = new double[_w * _h];
            _dst = new double[_w * _h];
        }

        private void swap() {
            _temp = _src;
            _src = _dst;
            _dst = _temp;
        }

        private double getSrc(int x, int y) {
            return _src[x + y * _w];
        }

        private void setSrc(int x, int y, double value) {
            _src[x + y * _w] = value;
        }

        /**
         * Linear interpolate on grid at coordinates (x, y). Coordinates will be
         * clamped to lie in simulation domain
         */
        private double lerpGrid(double x, double y) {
            x = min(max(x - _ox, 0.0), _w - 1.001);
            y = min(max(y - _oy, 0.0), _h - 1.001);
            int ix = (int) x;
            int iy = (int) y;
            x -= ix;
            y -= iy;

            double x00 = getSrc(ix + 0, iy + 0), x10 = getSrc(ix + 1, iy + 0);
            double x01 = getSrc(ix + 0, iy + 1), x11 = getSrc(ix + 1, iy + 1);

            return lerp(lerp(x00, x10, x), lerp(x01, x11, x), y);
        }

        /**
         * Advect grid in velocity field u, v with given timestep
         */
        private void advect(double timestep, FluidQuantity u, FluidQuantity v) {
            for (int iy = 0, idx = 0; iy < _h; iy++) {
                for (int ix = 0; ix < _w; ix++, idx++) {
                    double x = ix + _ox;
                    double y = iy + _oy;

                    /**
                     * First component: Integrate in time
                     */
                    Vector2D result = Euler(x, y, timestep, u, v);
                    x = result.getX();
                    y = result.getY();

                    /**
                     * Second component: Interpolate from grid
                     */
                    _dst[idx] = lerpGrid(x, y);
                }
            }
        }

        /**
         * Sets fluid quantity inside the given rect to value `v'
         */
        private void addInflow(double x0, double y0, double x1, double y1, double v) {
            int ix0 = (int) (x0 / _hx - _ox);
            int iy0 = (int) (y0 / _hx - _oy);
            int ix1 = (int) (x1 / _hx - _ox);
            int iy1 = (int) (y1 / _hx - _oy);

            for (int y = max(iy0, 0); y < min(iy1, _h); y++) {
                for (int x = max(ix0, 0); x < min(ix1, _h); x++) {
                    if (abs((float) _src[x + y * _w]) < abs((float) v)) {
                        _src[x + y * _w] = v;
                    }
                }
            }
        }

        /**
         * Simple forward Euler method for velocity integration in time.
         */
        private Vector2D Euler(double x, double y, double timestep, FluidQuantity u, FluidQuantity v) {
            double uVel = u.lerpGrid(x, y) / _hx;
            double vVel = v.lerpGrid(x, y) / _hx;

            x -= uVel * timestep;
            y -= vVel * timestep;

            return new Vector2D(x, y);
        }
    }

    /**
     * Fluid solver class. Sets up the fluid quantities, forces
     * incompressibility performs advection and adds inflows.
     */
    static class FluidSolver {

        /**
         * Fluid quantities
         */
        FluidQuantity _d;
        FluidQuantity _u;
        FluidQuantity _v;

        /**
         * Width and height
         */
        int _w;
        int _h;

        /**
         * Grid cell size and fluid density
         */
        double _hx;
        double _density;

        /**
         * Arrays for: Right hand side of pressure solve and pressure solution
         */
        double[] _r;
        double[] _p;

        FluidSolver(int w, int h, double density) {
            _w = w;
            _h = h;
            _density = density;
            _hx = 1.0 / min(w, h);

            _d = new FluidQuantity(_w, _h, 0.5, 0.5, _hx);
            _u = new FluidQuantity(_w + 1, _h, 0.0, 0.5, _hx);
            _v = new FluidQuantity(_w, _h + 1, 0.5, 0.0, _hx);

            _r = new double[_w * _h];
            _p = new double[_w * _h];
        }

        private void update(double timestep) {
            buildRhs();
            project(600, timestep);
            applyPressure(timestep);

            _d.advect(timestep, _u, _v);
            _u.advect(timestep, _u, _v);
            _v.advect(timestep, _u, _v);

            /**
             * Make effect of advection visible, since it's not an in-place
             * operation
             */
            _d.swap();
            _u.swap();
            _v.swap();
        }

        /**
         * Set density and x/y velocity in given rectangle to d/u/v,
         * respectively
         */
        private void addInflow(double x, double y, double w, double h, double d, double u, double v) {
            _d.addInflow(x, y, x + w, y + h, d);
            _u.addInflow(x, y, x + w, y + h, u);
            _v.addInflow(x, y, x + w, y + h, v);
        }

        /**
         * Returns the maximum allowed timestep. Note that the actual timestep
         * taken should usually be much below this to ensure accurate simulation
         * - just never above.
         */
        private double maxTimestep() {
            double maxVelocity = 0.0;
            for (int y = 0; y < _h; y++) {
                for (int x = 0; x < _w; x++) {
                    /**
                     * Average velocity at grid cell center
                     */
                    double u = _u.lerpGrid(x + 0.5, y + 0.5);
                    double v = _v.lerpGrid(x + 0.5, y + 0.5);

                    double velocity = sqrt(u * u + v * v);
                    maxVelocity = max(maxVelocity, velocity);
                }
            }

            /**
             * Fluid should not flow more than two grid cells per iteration
             */
            double maxTimestep = 2.0 * _hx / maxVelocity;

            /**
             * Clamp to sensible maximum value in case of very small velocities
             */
            return min(maxTimestep, 1.0);
        }

        /**
         * Convert fluid density to RGBA image.
         *
         * @param destination
         * @throws IOException
         */
        private void toImage(File destination) throws IOException {

            //write ppm header
            try (BufferedWriter output = new BufferedWriter(new FileWriter(destination))) {
                //write ppm header
                output.write("P3\n" + _w + " " + _h + "\n" + 255 + "\n");

                //for (int y = _w - 1; y >= 0; y--) {
                //    for (int x = 0; x < _h; x++) {
                for (int i = 0; i < _w * _h; i++) {
                    int shade = (int) ((1.0 - _d._src[i]) * 255.0);
                    shade = max(min(shade, 255), 0);
                    output.write(shade + " " + shade + " " + shade + " ");
                    if (i % _w == 0) {
                        output.write("\n");
                    }
                }
            }
        }

        /**
         * Builds the pressure right hand side as the negative divergence
         */
        private void buildRhs() {
            double scale = 1.0 / _hx;

            for (int y = 0, idx = 0; y < _h; y++) {
                for (int x = 0; x < _w; x++, idx++) {
                    _r[idx] = -scale * (_u.getSrc(x + 1, y) - _u.getSrc(x, y)
                            + _v.getSrc(x, y + 1) - _v.getSrc(x, y));
                }
            }
        }

        /**
         * Performs the pressure solve using Gauss-Seidel. The solver will run
         * as long as it takes to get the relative error below a threshold, but
         * will never exceed `limit' iterations
         */
        private void project(int limit, double timestep) {
            double scale = timestep / (_density * _hx * _hx);

            double maxDelta = 0.0;
            for (int iter = 0; iter < limit; iter++) {
                maxDelta = 0.0;
                for (int y = 0, idx = 0; y < _h; y++) {
                    for (int x = 0; x < _w; x++, idx++) {
                        double diag = 0.0, offDiag = 0.0;

                        /**
                         * Here we build the matrix implicitly as the five-point
                         * stencil. Grid borders are assumed to be solid, i.e.
                         * there is no fluid outside the simulation domain.
                         */
                        if (x > 0) {
                            diag += scale;
                            offDiag -= scale * _p[idx - 1];
                        }
                        if (y > 0) {
                            diag += scale;
                            offDiag -= scale * _p[idx - _w];
                        }
                        if (x < _w - 1) {
                            diag += scale;
                            offDiag -= scale * _p[idx + 1];
                        }
                        if (y < _h - 1) {
                            diag += scale;
                            offDiag -= scale * _p[idx + _w];
                        }

                        double newP = (_r[idx] - offDiag) / diag;

                        maxDelta = max(maxDelta, abs((float) _p[idx] - newP));

                        _p[idx] = newP;
                    }
                }

                if (maxDelta < 1e-5) {
                    log.log(INFO, format("Exiting solver after %d iterations, maximum change is %f", iter, maxDelta));
                    return;
                }
            }

            log.log(INFO, format("Exceeded budget of %d iterations, maximum change was %f", limit, maxDelta));
        }

        /**
         * Applies the computed pressure to the velocity field
         */
        private void applyPressure(double timestep) {
            double scale = timestep / (_density * _hx);

            for (int y = 0, idx = 0; y < _h; y++) {
                for (int x = 0; x < _w; x++, idx++) {
                    _u.setSrc(x, y, _u.getSrc(x, y) - scale * _p[idx]);
                    _u.setSrc(x + 1, y, _u.getSrc(x + 1, y) + scale * _p[idx]);
                    _v.setSrc(x, y, _v.getSrc(x, y) - scale * _p[idx]);
                    _v.setSrc(x, y + 1, _v.getSrc(x, y + 1) + scale * _p[idx]);
                }
            }

            for (int y = 0; y < _h; y++) {
                _u.setSrc(0, y, 0.0);
                _u.setSrc(_w, y, 0.0);
            }
            for (int x = 0; x < _w; x++) {
                _v.setSrc(x, 0, 0.0);
                _v.setSrc(x, _h, 0.0);
            }
        }
    }
}
