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
package physics.hybridFluids;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import static java.lang.Double.isNaN;
import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.random;
import static java.lang.Math.signum;
import static java.lang.Math.sqrt;
import static java.lang.String.format;
import static java.lang.System.arraycopy;
import static java.lang.System.exit;
import java.util.ArrayList;
import java.util.Stack;
import static java.util.logging.Level.ALL;
import static java.util.logging.Level.INFO;
import static java.util.logging.Level.SEVERE;
import java.util.logging.Logger;
import static java.util.logging.Logger.getLogger;
import static math.Interpolation.lerp;
import org.apache.commons.math3.geometry.euclidean.twod.Vector2D;
import static physics.hybridFluids.ExplicitFluid8flip.CellType.CELL_EMPTY;
import static physics.hybridFluids.ExplicitFluid8flip.CellType.CELL_FLUID;
import static physics.hybridFluids.ExplicitFluid8flip.CellType.CELL_SOLID;
import physics.solidBodies.A_SolidBody;
import physics.solidBodies.SolidBox;

/**
 * @author Gonçalo Amador {@literal &} Abel Gomes
 */
class ExplicitFluid8flip {

    static final Logger log = getLogger("ExplicitFluid8flip");

    static {
        log.setLevel(ALL);
    }

    public static void main(String args[]) {
        /**
         * Play with these constants, if you want
         */
        int sizeX = 128;
        int sizeY = 128;

        double densityAir = 0.1;
        double densitySoot = 0.25; /* You can make this smaller to get lighter smoke */

        double diffusion = 0.01;
        double timestep = 0.0025;

        boolean renderHeat = false; /* Set this to true to enable heat rendering */

        ArrayList<A_SolidBody> bodies = new ArrayList<>();
        bodies.add(new SolidBox(0.5, 0.6, 0.7, 0.1, PI * 0.25, 0.0, 0.0, 0.0));

        FluidSolver solver = new FluidSolver(sizeX, sizeY, densityAir, densitySoot, diffusion, bodies);

        double time = 0.0;
        int iterations = 0;

        while (time < 8.0) {
            for (int i = 0; i < 4; i++) {
                solver.update(timestep);
                time += timestep;
            }

            iterations++;
            File f = new File("Frame" + iterations + ".ppm");
            try {
                solver.toImage(f, renderHeat);
            } catch (IOException ex) {
                log.log(SEVERE, null, ex);
            }

            bodies.forEach(body -> body.update(timestep));
        }
    }

    /**
     * Length of ArrayList (x, y)
     *
     * @param x
     * @param y
     * @return
     */
    static private double length(double x, double y) {
        return sqrt(x * x + y * y);
    }

    /**
     * Cubic pulse function. Returns a value in range [0, 1]. Return value is 0
     * for x <= -1 and x >= 1; value is 1 for x=0 Smoothly interpolates between
     * 0 and 1 between these three points.
     *
     * @param x
     * @return
     */
    static private double cubicPulse(double x) {
        x = min(abs((float) x), 1.0);
        return 1.0 - x * x * (3.0 - 2.0 * x);
    }

    /**
     * For three corners in a 1x1 square, with 'in' being adjacent to 'out1' and
     * 'out2' and all three parameters being distances to a surface, 'in' being
     * inside the surface and 'out1' and 'out2' outside, returns the area of the
     * square occupied by the surface.
     */
    static private double triangleOccupancy(double out1, double in, double out2) {
        return 0.5 * in * in / ((out1 - in) * (out2 - in));
    }

    /**
     * For four corners in a 1x1 square, with all parameters being distances to
     * a surface and 'in1' and 'in2 inside the surface, returns the are of the
     * square occupied by the surface.
     */
    static private double trapezoidOccupancy(double out1, double out2, double in1, double in2) {
        return 0.5 * (-in1 / (out1 - in1) - in2 / (out2 - in2));
    }

    /**
     * Given the distance of four corners in a 1x1 square to a surface, returns
     * the area of the part of the square occupied by the surface computed
     * analytically.
     *
     * The basic workings of this algorithm are quite similar to marching
     * squares (2D marching cubes). First, a mask is computed based on which
     * corners are inside and which are outside. Based on this mask, the
     * function differentiates between one of four cases: a) Only one corner is
     * inside => Compute using triangle area b) Only one corner is outside =>
     * Invert distance field, compute 1 - triangle area c) Two adjacent corners
     * are inside => Compute using trapezoid area d) Two opposing corners are
     * inside => Compute as sum of area of two opposed triangles
     *
     * The two remaining cases, all corners outside/inside, can be computed
     * trivially
     */
    static private double occupancy(double d11, double d12, double d21, double d22) {
        double ds[] = {d11, d12, d22, d21};

        /**
         * Compute mask
         */
        int b = 0;
        for (int i = 3; i >= 0; i--) {
            b = (b << 1) | (ds[i] < 0.0 ? 1 : 0);
        }

        switch (b) {
            case 0x0:
                return 0.0;

            case 0x1:
                return triangleOccupancy(d21, d11, d12);
            case 0x2:
                return triangleOccupancy(d11, d12, d22);
            case 0x4:
                return triangleOccupancy(d12, d22, d21);
            case 0x8:
                return triangleOccupancy(d22, d21, d11);

            case 0xE:
                return 1.0 - triangleOccupancy(-d21, -d11, -d12);
            case 0xD:
                return 1.0 - triangleOccupancy(-d11, -d12, -d22);
            case 0xB:
                return 1.0 - triangleOccupancy(-d12, -d22, -d21);
            case 0x7:
                return 1.0 - triangleOccupancy(-d22, -d21, -d11);

            case 0x3:
                return trapezoidOccupancy(d21, d22, d11, d12);
            case 0x6:
                return trapezoidOccupancy(d11, d21, d12, d22);
            case 0x9:
                return trapezoidOccupancy(d12, d22, d11, d21);
            case 0xC:
                return trapezoidOccupancy(d11, d12, d21, d22);

            case 0x5:
                return triangleOccupancy(d11, d12, d22)
                        + triangleOccupancy(d22, d21, d11);
            case 0xA:
                return triangleOccupancy(d21, d11, d12)
                        + triangleOccupancy(d12, d22, d21);

            case 0xF:
                return 1.0;
        }

        return 0.0;
    }

    static class FluidQuantity {

        double[] _src;
        double[] _old; /* Contains old quantities at beginning of iteration */

        /**
         * Distance field induced by solids. Since this is used to compute the
         * cell volumes, the samples are offset by (-0.5, -0.5) from the samples
         * in _src and the grid is one larger in each dimension. This way, each
         * sample of fluid quantity has four samples of the distance function
         * surrounding it - perfect for computing the cell volumes.
         */
        double[] _phi;

        /**
         * Fractional cell volume occupied by fluid
         */
        double[] _volume;

        /**
         * Normal of distance field at grid points
         */
        double[] _normalX;
        double[] _normalY;

        /**
         * Designates cells as fluid or solid cells (CELL_FLUID or CELL_SOLID)
         */
        CellType[] _cell;

        /**
         * Specifies the index of the rigid body closes to a grid cell
         */
        int[] _body;
        /**
         * Auxiliary array used for extrapolation
         */
        int[] _mask;

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
        double _hx;

        FluidQuantity(int w, int h, double ox, double oy, double hx) {
            _w = w;
            _h = h;
            _ox = ox;
            _oy = oy;
            _hx = hx;

            _src = new double[_w * _h];
            _old = new double[_w * _h];

            /**
             * Make distance grid one larger in each dimension
             */
            _phi = new double[(_w + 1) * (_h + 1)];
            _volume = new double[_w * _h];
            _normalX = new double[_w * _h];
            _normalY = new double[_w * _h];

            _cell = new CellType[_w * _h];

            _body = new int[_w * _h];
            _mask = new int[_w * _h];

            for (int i = 0; i < _w * _h; i++) {
                _cell[i] = CELL_FLUID;
                _volume[i] = 1.0;
            }
        }

        /**
         * Adds contribution 'value' of sample at (x, y) to grid cell at (ix,
         * iy) using a hat filter.
         */
        private void addSample(double[] weight, double value, double x, double y, int ix, int iy) {
            if (ix < 0 || iy < 0 || ix >= _w || iy >= _h) {
                return;
            }

            //double k = (1.0 - abs((float) ix - x)) * (1.0 - abs((float) iy - y));
            double k = (1.0 - abs(ix - x)) * (1.0 - abs(iy - y));
            weight[ix + iy * _w] += k;
            _src[ix + iy * _w] += k * value;
        }

        private int I(int x, int y) {
            return x + y * _w;
        }

        private double getSrc(int x, int y) {
            return _src[x + y * _w];
        }

        private double volume(int x, int y) {
            return _volume[x + y * _w];
        }

        private void setSrc(int x, int y, double value) {
            _src[x + y * _w] = value;
        }

        private void copy() {
            arraycopy(_src, 0, _old, 0, _w * _h);
        }

        /**
         * Computes the change in quantity during the last update
         */
        private void diff(double alpha) {
            for (int i = 0; i < _w * _h; i++) {
                _src[i] -= (1.0 - alpha) * _old[i];
            }
        }

        /**
         * Reverses the previous transformation - saves memory
         */
        private void undiff(double alpha) {
            for (int i = 0; i < _w * _h; i++) {
                _src[i] += (1.0 - alpha) * _old[i];
            }
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
         * Set fluid quantity inside the given rect to the specified value, but
         * use a smooth falloff to avoid oscillations
         */
        private void addInflow(double x0, double y0, double x1, double y1, double v) {
            int ix0 = (int) (x0 / _hx - _ox);
            int iy0 = (int) (y0 / _hx - _oy);
            int ix1 = (int) (x1 / _hx - _ox);
            int iy1 = (int) (y1 / _hx - _oy);

            for (int y = max(iy0, 0); y < min(iy1, _h); y++) {
                for (int x = max(ix0, 0); x < min(ix1, _h); x++) {
                    double l = length(
                            (2.0 * (x + 0.5) * _hx - (x0 + x1)) / (x1 - x0),
                            (2.0 * (y + 0.5) * _hx - (y0 + y1)) / (y1 - y0)
                    );
                    double vi = cubicPulse(l) * v;
                    if (abs((float) _src[x + y * _w]) < abs((float) vi)) {
                        _src[x + y * _w] = vi;
                    }
                }
            }
        }

        /**
         * Fill all solid related fields - that is, _cell, _body and _normalX/Y
         */
        private void fillSolidFields(ArrayList<A_SolidBody> bodies) {
            if (bodies.isEmpty()) {
                return;
            }

            /**
             * Compute distance field first
             */
            for (int iy = 0, idx = 0; iy <= _h; iy++) {
                for (int ix = 0; ix <= _w; ix++, idx++) {
                    double x = (ix + _ox - 0.5) * _hx;
                    double y = (iy + _oy - 0.5) * _hx;

                    _phi[idx] = bodies.get(0).distance(x, y);
                    for (int i = 1; i < bodies.size(); i++) {
                        _phi[idx] = min(_phi[idx], bodies.get(i).distance(x, y));
                    }
                }
            }

            for (int iy = 0, idx = 0; iy < _h; iy++) {
                for (int ix = 0; ix < _w; ix++, idx++) {
                    double x = (ix + _ox) * _hx;
                    double y = (iy + _oy) * _hx;

                    _body[idx] = 0;
                    double d = bodies.get(0).distance(x, y);
                    for (int i = 1; i < bodies.size(); i++) {
                        double id = bodies.get(i).distance(x, y);
                        if (id < d) {
                            _body[idx] = i;
                            d = id;
                        }
                    }

                    /**
                     * Compute cell volume from the four adjacent distance
                     * samples
                     */
                    int idxp = ix + iy * (_w + 1);
                    _volume[idx] = 1.0 - occupancy(
                            _phi[idxp], _phi[idxp + 1],
                            _phi[idxp + _w + 1], _phi[idxp + _w + 2]
                    );

                    /**
                     * Clamp dangerously small cell volumes - could break
                     * numerical solver otherwise
                     */
                    if (_volume[idx] < 0.01) {
                        _volume[idx] = 0.0;
                    }

                    /**
                     * Solid cells are now defined as cells with zero fluid
                     * volume
                     */
                    Vector2D result = bodies.get(_body[idx]).distanceNormal(_normalX[idx], _normalY[idx], x, y);
                    _normalX[idx] = result.getX();
                    _normalY[idx] = result.getY();

                    if (_volume[idx] == 0.0) {
                        _cell[idx] = CELL_SOLID;
                    } else {
                        _cell[idx] = CELL_FLUID;
                    }
                }
            }
        }

        /**
         * The extrapolation routine is now augmented to also fill in values for
         * cells that ended up with no particles in them. These are marked with
         * CELL_EMPTY. Empty cells are computed as the average value of all
         * available neighbours, and can therefore be computed as soon as at
         * least one neighbouring cell is available.
         */
        private void fillSolidMask() {
            /**
             * Make sure border is not touched by extrapolation - will be
             * handled separately.
             */
            for (int x = 0; x < _w; x++) {
                _mask[x] = _mask[x + (_h - 1) * _w] = 0xFF;
            }
            for (int y = 0; y < _h; y++) {
                _mask[y * _w] = _mask[y * _w + _w - 1] = 0xFF;
            }

            for (int y = 1; y < _h - 1; y++) {
                for (int x = 1; x < _w - 1; x++) {
                    int idx = x + y * _w;

                    _mask[idx] = 0;
                    if (_cell[idx] == CELL_SOLID) {
                        double nx = _normalX[idx];
                        double ny = _normalY[idx];

                        if (nx != 0.0 && _cell[idx + (int) signum(nx)] != CELL_FLUID) {
                            _mask[idx] |= 1;
                        }
                        if (ny != 0.0 && _cell[idx + (int) signum(ny) * _w] != CELL_FLUID) {
                            _mask[idx] |= 2;
                        }
                    } else if (_cell[idx] == CELL_EMPTY) {
                        /**
                         * Empty cells with no available neighbours need to be
                         * processed later.
                         */
                        _mask[idx]
                                = _cell[idx - 1] != CELL_FLUID
                                && _cell[idx + 1] != CELL_FLUID
                                && _cell[idx - _w] != CELL_FLUID
                                && _cell[idx + _w] != CELL_FLUID ? 1 : 0;
                    }
                }
            }
        }

        /**
         * Solve for value at index idx using values of neighbours in normal x/y
         * direction. The value is computed such that the directional derivative
         * along distance field normal is 0.
         */
        private double extrapolateNormal(int idx) {
            double nx = _normalX[idx];
            double ny = _normalY[idx];

            double srcX = _src[idx + (int) signum(nx)];
            double srcY = _src[idx + (int) signum(ny) * _w];

            return (abs((float) nx) * srcX + abs((float) ny) * srcY) / (abs((float) nx) + abs((float) ny));
        }

        /**
         * Computes the extrapolated value as the average of all available
         * neighbouring cells.
         */
        private double extrapolateAverage(int idx) {
            double value = 0.0;
            int count = 0;

            if (_cell[idx - 1] == CELL_FLUID) {
                value += _src[idx - 1];
                count++;
            }
            if (_cell[idx + 1] == CELL_FLUID) {
                value += _src[idx + 1];
                count++;
            }
            if (_cell[idx - _w] == CELL_FLUID) {
                value += _src[idx - _w];
                count++;
            }
            if (_cell[idx + _w] == CELL_FLUID) {
                value += _src[idx + _w];
                count++;
            }
            return value / count;
        }

        private Stack<Integer> freeSolidNeighbour(
                int idx,
                Stack<Integer> border,
                int mask) {
            if (_cell[idx] == CELL_SOLID) {
                _mask[idx] &= ~mask;
                if (_mask[idx] == 0) {
                    border.push(idx);
                }
            }
            return border;
        }

        /**
         * At least one free neighbour cell is enough to add this cell to the
         * queue of ready cells.
         */
        private Stack<Integer> freeEmptyNeighbour(
                int idx,
                Stack<Integer> border) {
            if (_cell[idx] == CELL_EMPTY && _mask[idx] == 1) {
                _mask[idx] = 0;
                border.push(idx);
            }
            return border;
        }

        /**
         * For empty cells on the border of the simulation domain, we simply
         * copy the values of the adjacent cells.
         */
        private void extrapolateEmptyBorders() {
            for (int x = 1; x < _w - 1; x++) {
                int idxT = x;
                int idxB = x + (_h - 1) * _w;

                if (_cell[idxT] == CELL_EMPTY) {
                    _src[idxT] = _src[idxT + _w];
                }
                if (_cell[idxB] == CELL_EMPTY) {
                    _src[idxB] = _src[idxB - _w];
                }
            }

            for (int y = 1; y < _h - 1; y++) {
                int idxL = y * _w;
                int idxR = y * _w + _w - 1;

                if (_cell[idxL] == CELL_EMPTY) {
                    _src[idxL] = _src[idxL + 1];
                }
                if (_cell[idxR] == CELL_EMPTY) {
                    _src[idxR] = _src[idxR - 1];
                }
            }

            int idxTL = 0;
            int idxTR = _w - 1;
            int idxBL = (_h - 1) * _w;
            int idxBR = _h * _w - 1;

            /**
             * Corner cells average the values of the two adjacent border cells
             */
            if (_cell[idxTL] == CELL_EMPTY) {
                _src[idxTL] = 0.5 * (_src[idxTL + 1] + _src[idxTL + _w]);
            }
            if (_cell[idxTR] == CELL_EMPTY) {
                _src[idxTR] = 0.5 * (_src[idxTR - 1] + _src[idxTR + _w]);
            }
            if (_cell[idxBL] == CELL_EMPTY) {
                _src[idxBL] = 0.5 * (_src[idxBL + 1] + _src[idxBL - _w]);
            }
            if (_cell[idxBR] == CELL_EMPTY) {
                _src[idxBR] = 0.5 * (_src[idxBR - 1] + _src[idxBR - _w]);
            }

            for (int i = 0; i < _w * _h; i++) {
                if (_cell[i] == CELL_EMPTY) {
                    _cell[i] = CELL_FLUID;
                }
            }
        }

        private void extrapolate() {
            fillSolidMask();

            /**
             * Queue of cells which can be computed
             */
            Stack<Integer> border = new Stack<>();

            /**
             * Initialize queue by finding all solid cells with mask=0 (ready
             * for extrapolation)
             */
            for (int y = 1; y < _h - 1; y++) {
                for (int x = 1; x < _w - 1; x++) {
                    int idx = x + y * _w;

                    if (_cell[idx] != CELL_FLUID && _mask[idx] == 0) {
                        border.push(idx);
                    }
                }
            }

            while (!border.empty()) {
                int idx = border.peek();
                border.pop();

                if (_cell[idx] == CELL_EMPTY) {
                    _src[idx] = extrapolateAverage(idx);
                    _cell[idx] = CELL_FLUID; /* Mark extrapolated empty cells as fluid */

                } else {
                    _src[idx] = extrapolateNormal(idx);
                }

                if (_normalX[idx - 1] > 0.0) {
                    border = freeSolidNeighbour(idx - 1, border, 1);
                }
                if (_normalX[idx + 1] < 0.0) {
                    border = freeSolidNeighbour(idx + 1, border, 1);
                }
                if (_normalY[idx - _w] > 0.0) {
                    border = freeSolidNeighbour(idx - _w, border, 2);
                }
                if (_normalY[idx + _w] < 0.0) {
                    border = freeSolidNeighbour(idx + _w, border, 2);
                }

                /**
                 * Notify adjacent empty cells
                 */
                border = freeEmptyNeighbour(idx - 1, border);
                border = freeEmptyNeighbour(idx + 1, border);
                border = freeEmptyNeighbour(idx - _w, border);
                border = freeEmptyNeighbour(idx + _w, border);
            }

            extrapolateEmptyBorders();
        }

        /**
         * Transfers particle values onto grid using a linear filter.
         *
         * In a first step, particle values and filter weights are accumulated
         * on the grid by looping over all particles and adding the particle
         * contribution to the four closest grid cells.
         *
         * In a second step, the actual grid values are obtained by dividing by
         * the filter weights. Cells with weight zero are cells which do not
         * contain any particles and are subsequently marked as empty for
         * extrapolation.
         */
        private void fromParticles(
                double[] weight,
                int count,
                double[] posX,
                double[] posY,
                double[] property) {
            _src = new double[_w * _h];
            weight = new double[_w * _h];

            for (int i = 0; i < count; i++) {
                double x = posX[i] - _ox;
                double y = posY[i] - _oy;
                x = max(0.5, min(_w - 1.5, x));
                y = max(0.5, min(_h - 1.5, y));

                int ix = (int) x;
                int iy = (int) y;

                addSample(weight, property[i], x, y, ix + 0, iy + 0);
                addSample(weight, property[i], x, y, ix + 1, iy + 0);
                addSample(weight, property[i], x, y, ix + 0, iy + 1);
                addSample(weight, property[i], x, y, ix + 1, iy + 1);
            }

            for (int i = 0; i < _w * _h; i++) {
                if (weight[i] != 0.0) {
                    _src[i] /= weight[i];
                } else if (_cell[i] == CELL_FLUID) {
                    _cell[i] = CELL_EMPTY;
                }
            }
        }
    }

    /**
     * Fluid solver class. Sets up the fluid quantities, forces
     * incompressibility performs advection and adds inflows.
     */
    static class FluidSolver {

        static ParticleQuantities _qs;

        /**
         * Width and height
         */
        static int _w;
        static int _h;

        /**
         * Grid cell size
         */
        static double _hx;

        static ArrayList<A_SolidBody> _bodies;

        /**
         * Fluid quantities
         */
        static FluidQuantity _d;
        static FluidQuantity _t; /* Temperature */

        static FluidQuantity _u;
        static FluidQuantity _v;

        /**
         * Densities at staggered grid locations
         */
        double[] _uDensity;
        double[] _vDensity;

        /**
         * Density of air
         */
        double _densityAir;

        /**
         * Density of soot
         */
        double _densitySoot;

        /**
         * Diffusion rate of heat
         */
        double _diffusion;

        /**
         * Arrays for: Right hand side of pressure solve and pressure solution
         */
        double[] _r;
        double[] _p;

        /**
         * Auxiliary vector
         */
        double[] _z;

        /**
         * Search vector
         */
        double[] _s;

        /**
         * Preconditioner
         */
        double[] _precon;

        /**
         * Matrix diagonal
         */
        double[] _aDiag;

        /**
         * Matrix off-X-diagonal
         */
        double[] _aPlusX;

        /**
         * Matrix off-Y-diagonal
         */
        double[] _aPlusY;

        /**
         * Ambient temperature (here room temperature), in Kelvin
         */
        double _tAmb;

        /**
         * Gravity
         */
        double _g;

        /**
         * Tiny blending factor for FLIP/PIC to avoid noise
         */
        double _flipAlpha;

        FluidSolver(int w, int h, double rhoAir, double rhoSoot, double diffusion,
                ArrayList<A_SolidBody> bodies) {
            _w = w;
            _h = h;
            _densityAir = rhoAir;
            _densitySoot = rhoSoot;
            _diffusion = diffusion;
            _bodies = bodies;

            _tAmb = 294.0;
            _g = 9.81;
            _flipAlpha = 0.001;

            _hx = 1.0 / min(w, h);

            _d = new FluidQuantity(_w, _h, 0.5, 0.5, _hx);
            _t = new FluidQuantity(_w, _h, 0.5, 0.5, _hx);
            _u = new FluidQuantity(_w + 1, _h, 0.0, 0.5, _hx);
            _v = new FluidQuantity(_w, _h + 1, 0.5, 0.0, _hx);

            for (int i = 0; i < _w * _h; i++) {
                _t._src[i] = _tAmb;
            }

            _qs = new ParticleQuantities(_w, _h, _hx, _bodies);
            _qs.addQuantity(_d);
            _qs.addQuantity(_t);
            _qs.addQuantity(_u);
            _qs.addQuantity(_v);

            /**
             * Interpolate initial quantity distribution onto particles
             */
            _qs.gridToParticles(1.0);

            _r = new double[_w * _h];
            _p = new double[_w * _h];
            _z = new double[_w * _h];
            _s = new double[_w * _h];
            _aDiag = new double[_w * _h];
            _aPlusX = new double[_w * _h];
            _aPlusY = new double[_w * _h];
            _precon = new double[_w * _h];

            _uDensity = new double[(_w + 1) * _h];
            _vDensity = new double[_w * (_h + 1)];
        }

        /**
         * We now modify the right hand side to "blend" between solid and fluid
         * velocity based on the cell volume occupied by fluid.
         */
        private void buildRhs() {
            double scale = 1.0 / _hx;
            CellType[] cell = _d._cell;
            int[] body = _d._body;

            for (int y = 0, idx = 0; y < _h; y++) {
                for (int x = 0; x < _w; x++, idx++) {
                    if (cell[idx] == CELL_FLUID) {
                        _r[idx] = -scale
                                * (_u.volume(x + 1, y) * _u.getSrc(x + 1, y) - _u.volume(x, y) * _u.getSrc(x, y)
                                + _v.volume(x, y + 1) * _v.getSrc(x, y + 1) - _v.volume(x, y) * _v.getSrc(x, y));

                        double vol = _d.volume(x, y);

                        if (_bodies.isEmpty()) {
                            continue;
                        }

                        if (x > 0) {
                            _r[idx] -= (_u.volume(x, y) - vol) * _bodies.get(body[idx - 1]).velocityX(x * _hx, (y + 0.5) * _hx);
                        }
                        if (y > 0) {
                            _r[idx] -= (_v.volume(x, y) - vol) * _bodies.get(body[idx - _w]).velocityY((x + 0.5) * _hx, y * _hx);
                        }
                        if (x < _w - 1) {
                            _r[idx] += (_u.volume(x + 1, y) - vol) * _bodies.get(body[idx + 1]).velocityX((x + 1.0) * _hx, (y + 0.5) * _hx);
                        }
                        if (y < _h - 1) {
                            _r[idx] += (_v.volume(x, y + 1) - vol) * _bodies.get(body[idx + _w]).velocityY((x + 0.5) * _hx, (y + 1.0) * _hx);
                        }
                    } else {
                        _r[idx] = 0.0;
                    }
                }
            }
        }

        /**
         * Computes densities at the staggered grid locations as a function of
         * temperature and smoke concentration.
         */
        private void computeDensities() {
            double alpha = (_densitySoot - _densityAir) / _densityAir;

            _uDensity = new double[(_w + 1) * _h];
            _vDensity = new double[_w * (_h + 1)];

            for (int y = 0; y < _h; y++) {
                for (int x = 0; x < _w; x++) {
                    double density = _densityAir * _tAmb / _t.getSrc(x, y) * (1.0 + alpha * _d.getSrc(x, y));
                    density = max(density, 0.05 * _densityAir);

                    _uDensity[_u.I(x, y)] += 0.5 * density;
                    _vDensity[_v.I(x, y)] += 0.5 * density;
                    _uDensity[_u.I(x + 1, y)] += 0.5 * density;
                    _vDensity[_v.I(x, y + 1)] += 0.5 * density;
                }
            }
        }

        /**
         * Instead of ant density per cell, the entries must now be modified to
         * account for variable density at individual grid cells.
         */
        private void buildPressureMatrix(double timestep) {
            double scale = timestep / (_hx * _hx);
            CellType[] cell = _d._cell;

            _aDiag = new double[_w * _h];
            _aPlusX = new double[_w * _h];
            _aPlusY = new double[_w * _h];

            for (int y = 0, idx = 0; y < _h; y++) {
                for (int x = 0; x < _w; x++, idx++) {
                    if (cell[idx] != CELL_FLUID) {
                        continue;
                    }

                    if (x < _w - 1 && cell[idx + 1] == CELL_FLUID) {
                        double factor = scale * _u.volume(x + 1, y) / _uDensity[_u.I(x + 1, y)];
                        _aDiag[idx] += factor;
                        _aDiag[idx + 1] += factor;
                        _aPlusX[idx] = -factor;
                    }
                    if (y < _h - 1 && cell[idx + _w] == CELL_FLUID) {
                        double factor = scale * _v.volume(x, y + 1) / _vDensity[_u.I(x, y + 1)];
                        _aDiag[idx] += factor;
                        _aDiag[idx + _w] += factor;
                        _aPlusY[idx] = -factor;
                    }
                }
            }
        }

        /**
         * Constructs the matrix used to compute heat diffusion. This uses the
         * diffusion ant, which specifies how quickly heat diffuses in the
         * domain. Higher values imply faster diffusion and a more difficult
         * systems of equations (CG will converge slower)
         */
        private void buildHeatDiffusionMatrix(double timestep) {
            for (int i = 0; i < _w * _h; i++) {
                _aDiag[i] = 1.0;
            }

            _aPlusX = new double[_w * _h];
            _aPlusY = new double[_w * _h];

            CellType[] cell = _d._cell;
            double scale = _diffusion * timestep * 1.0 / (_hx * _hx);

            for (int y = 0, idx = 0; y < _h; y++) {
                for (int x = 0; x < _w; x++, idx++) {
                    if (cell[idx] != CELL_FLUID) {
                        continue;
                    }

                    if (x < _w - 1 && cell[idx + 1] == CELL_FLUID) {
                        _aDiag[idx] += scale;
                        _aDiag[idx + 1] += scale;
                        _aPlusX[idx] = -scale;
                    }

                    if (y < _h - 1 && cell[idx + _w] == CELL_FLUID) {
                        _aDiag[idx] += scale;
                        _aDiag[idx + _w] += scale;
                        _aPlusY[idx] = -scale;
                    }
                }
            }
        }

        /**
         * Builds the modified incomplete Cholesky preconditioner
         */
        private void buildPreconditioner() {
            double tau = 0.97;
            double sigma = 0.25;
            CellType[] cell = _d._cell;

            for (int y = 0, idx = 0; y < _h; y++) {
                for (int x = 0; x < _w; x++, idx++) {
                    if (cell[idx] != CELL_FLUID) {
                        continue;
                    }

                    double e = _aDiag[idx];

                    if (x > 0 && cell[idx - 1] == CELL_FLUID) {
                        double px = _aPlusX[idx - 1] * _precon[idx - 1];
                        double py = _aPlusY[idx - 1] * _precon[idx - 1];
                        e -= (px * px + tau * px * py);
                    }
                    if (y > 0 && cell[idx - _w] == CELL_FLUID) {
                        double px = _aPlusX[idx - _w] * _precon[idx - _w];
                        double py = _aPlusY[idx - _w] * _precon[idx - _w];
                        e -= (py * py + tau * px * py);
                    }

                    if (e < sigma * _aDiag[idx]) {
                        e = _aDiag[idx];
                    }

                    _precon[idx] = 1.0 / sqrt(e);
                }
            }
        }

        /**
         * Apply preconditioner to vector 'a' and store it in 'dst'
         */
        private void applyPreconditioner(double[] dst, double[] a) {
            CellType[] cell = _d._cell;

            for (int y = 0, idx = 0; y < _h; y++) {
                for (int x = 0; x < _w; x++, idx++) {
                    if (cell[idx] != CELL_FLUID) {
                        continue;
                    }

                    double t = a[idx];

                    if (x > 0 && cell[idx - 1] == CELL_FLUID) {
                        t -= _aPlusX[idx - 1] * _precon[idx - 1] * dst[idx - 1];
                    }
                    if (y > 0 && cell[idx - _w] == CELL_FLUID) {
                        t -= _aPlusY[idx - _w] * _precon[idx - _w] * dst[idx - _w];
                    }

                    dst[idx] = t * _precon[idx];
                }
            }

            for (int y = _h - 1, idx = _w * _h - 1; y >= 0; y--) {
                for (int x = _w - 1; x >= 0; x--, idx--) {
                    if (cell[idx] != CELL_FLUID) {
                        continue;
                    }

                    double t = dst[idx];

                    if (x < _w - 1 && cell[idx + 1] == CELL_FLUID) {
                        t -= _aPlusX[idx] * _precon[idx] * dst[idx + 1];
                    }
                    if (y < _h - 1 && cell[idx + _w] == CELL_FLUID) {
                        t -= _aPlusY[idx] * _precon[idx] * dst[idx + _w];
                    }

                    dst[idx] = t * _precon[idx];
                }
            }
        }

        /**
         * Returns the dot product of vectors 'a' and 'b'
         */
        private double dotProduct(double[] a, double[] b) {
            CellType[] cell = _d._cell;

            double result = 0.0;
            for (int i = 0; i < _w * _h; i++) {
                if (cell[i] == CELL_FLUID) {
                    result += a[i] * b[i];
                }
            }
            return result;
        }

        /**
         * Multiplies internal pressure matrix with vector 'b' and stores the
         * result in 'dst'
         */
        private void matrixVectorProduct(double[] dst, double[] b) {
            for (int y = 0, idx = 0; y < _h; y++) {
                for (int x = 0; x < _w; x++, idx++) {
                    double t = _aDiag[idx] * b[idx];

                    if (x > 0) {
                        t += _aPlusX[idx - 1] * b[idx - 1];
                    }
                    if (y > 0) {
                        t += _aPlusY[idx - _w] * b[idx - _w];
                    }
                    if (x < _w - 1) {
                        t += _aPlusX[idx] * b[idx + 1];
                    }
                    if (y < _h - 1) {
                        t += _aPlusY[idx] * b[idx + _w];
                    }

                    dst[idx] = t;
                }
            }
        }

        /**
         * Computes 'dst' = 'a' + 'b'*'s'
         */
        private void scaledAdd(double[] dst, double[] a, double[] b, double s) {
            CellType[] cell = _d._cell;

            for (int i = 0; i < _w * _h; i++) {
                if (cell[i] == CELL_FLUID) {
                    dst[i] = a[i] + b[i] * s;
                }
            }
        }

        /**
         * Returns maximum absolute value in vector 'a'
         */
        private double infinityNorm(double[] a) {
            CellType[] cell = _d._cell;

            double maxA = 0.0;
            for (int i = 0; i < _w * _h; i++) {
                if (cell[i] == CELL_FLUID) {
                    maxA = max(maxA, abs((float) a[i]));
                }
            }
            return maxA;
        }

        /**
         * Conjugate gradients solver
         */
        private void project(int limit) {
            /**
             * Initial guess of zeroes
             */
            _p = new double[_w * _h];
            applyPreconditioner(_z, _r);
            arraycopy(_z, 0, _s, 0, _w * _h);

            double maxError = infinityNorm(_r);
            if (maxError < 1e-5) {
                log.info("Initial guess sufficiently small");
                return;
            }

            double sigma = dotProduct(_z, _r);

            for (int iter = 0; iter < limit; iter++) {
                matrixVectorProduct(_z, _s);
                double alpha = sigma / dotProduct(_z, _s);
                scaledAdd(_p, _p, _s, alpha);
                scaledAdd(_r, _r, _z, -alpha);

                maxError = infinityNorm(_r);
                if (maxError < 1e-5) {
                    log.log(INFO, format("Exiting solver after %d iterations, maximum error is %f", iter, maxError));
                    return;
                }
                if (isNaN(maxError)) {
                    log.log(INFO, format("Exiting solver after %d iterations, maximum error is %f", iter, maxError));
                    exit(-1);
                }

                applyPreconditioner(_z, _r);

                double sigmaNew = dotProduct(_z, _r);
                scaledAdd(_s, _z, _s, sigmaNew / sigma);
                sigma = sigmaNew;
            }

            log.log(INFO, format("Exceeded budget of %d iterations, maximum error was %f", limit, maxError));
        }

        /**
         * Similar to the pressure matrix, we cannot assume ant density per cell
         * here either and must modify the equations accordingly.
         */
        private void applyPressure(double timestep) {
            double scale = timestep / _hx;
            CellType[] cell = _d._cell;

            for (int y = 0, idx = 0; y < _h; y++) {
                for (int x = 0; x < _w; x++, idx++) {
                    if (cell[idx] != CELL_FLUID) {
                        continue;
                    }

                    _u.setSrc(x, y, _u.getSrc(x, y) - scale * _p[idx] / _uDensity[_u.I(x, y)]);
                    _v.setSrc(x, y, _v.getSrc(x, y) - scale * _p[idx] / _vDensity[_v.I(x, y)]);
                    _u.setSrc(x + 1, y, _u.getSrc(x + 1, y) + scale * _p[idx] / _uDensity[_u.I(x + 1, y)]);
                    _v.setSrc(x, y + 1, _v.getSrc(x, y + 1) + scale * _p[idx] / _vDensity[_v.I(x, y + 1)]);
                }
            }
        }

        /**
         * Add body force due to density and heat difference
         */
        private void addBuoyancy(double timestep) {
            double alpha = (_densitySoot - _densityAir) / _densityAir;

            for (int y = 0; y < _h; y++) {
                for (int x = 0; x < _w; x++) {
                    double buoyancy = timestep * _g * (alpha * _d.getSrc(x, y) - (_t.getSrc(x, y) - _tAmb) / _tAmb);

                    _v.setSrc(x, y, _v.getSrc(x, y) + buoyancy * 0.5);
                    _v.setSrc(x, y + 1, _v.getSrc(x, y + 1) + buoyancy * 0.5);
                }
            }
        }

        /**
         * Sets all velocity cells bordering solid cells to the solid velocity
         */
        private void setBoundaryCondition() {
            CellType[] cell = _d._cell;
            int[] body = _d._body;

            for (int y = 0, idx = 0; y < _h; y++) {
                for (int x = 0; x < _w; x++, idx++) {
                    if (cell[idx] == CELL_SOLID) {
                        A_SolidBody b = _bodies.get(body[idx]);

                        _u.setSrc(x, y, b.velocityX(x * _hx, (y + 0.5) * _hx));
                        _v.setSrc(x, y, b.velocityX((x + 0.5) * _hx, y * _hx));
                        _u.setSrc(x + 1, y, b.velocityX((x + 1.0) * _hx, (y + 0.5) * _hx));
                        _v.setSrc(x, y + 1, b.velocityX((x + 0.5) * _hx, (y + 1.0) * _hx));
                    }
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

        private void update(double timestep) {
            _d.fillSolidFields(_bodies);
            _t.fillSolidFields(_bodies);
            _u.fillSolidFields(_bodies);
            _v.fillSolidFields(_bodies);

            /**
             * Interpolate particle quantities to grid
             */
            _qs.particlesToGrid();

            /**
             * Set current values as the old/pre-update values
             */
            _d.copy();
            _t.copy();
            _u.copy();
            _v.copy();

            /**
             * Unfortunately, we have to move inflows out of the mainloop into
             * here - all changes need to happen between copy and diff to have
             * any effect
             */
            addInflow(0.45, 0.2, 0.2, 0.05, 1.0, _tAmb, 0.0, 0.0);

            /**
             * Right-hand side of heat equation is the current heat distribution
             */
            arraycopy(_t._src, 0, _r, 0, _w * _h);
            buildHeatDiffusionMatrix(timestep);
            buildPreconditioner();
            project(2_000);

            /**
             * The solution of the heat equation is the heat distribution in the
             * next timestep.
             */
            arraycopy(_p, 0, _t._src, 0, _w * _h);

            _t.extrapolate();

            addBuoyancy(timestep);
            setBoundaryCondition();

            buildRhs();
            computeDensities();
            buildPressureMatrix(timestep);
            buildPreconditioner();
            project(2_000);
            applyPressure(timestep);

            _d.extrapolate();
            _u.extrapolate();
            _v.extrapolate();

            setBoundaryCondition();

            /**
             * Compute change in quantities
             */
            _d.diff(_flipAlpha);
            _t.diff(_flipAlpha);
            _u.diff(_flipAlpha);
            _v.diff(_flipAlpha);

            /**
             * Interpolate change onto particles
             */
            _qs.gridToParticles(_flipAlpha);

            /**
             * Reverse the change computation to get the post-update values back
             * (for rendering/advection).
             */
            _d.undiff(_flipAlpha);
            _t.undiff(_flipAlpha);
            _u.undiff(_flipAlpha);
            _v.undiff(_flipAlpha);

            /**
             * Advect particles in velocity field
             */
            _qs.advect(timestep, _u, _v);
        }

        /**
         * Set density and x/y velocity in given rectangle to d/u/v,
         * respectively
         */
        private void addInflow(
                double x, double y,
                double w, double h,
                double d, double t,
                double u, double v) {
            _d.addInflow(x, y, x + w, y + h, d);
            _t.addInflow(x, y, x + w, y + h, t);
            _u.addInflow(x, y, x + w, y + h, u);
            _v.addInflow(x, y, x + w, y + h, v);
        }

        private double ambientT() {
            return _tAmb;
        }

        /**
         * Convert fluid density to RGBA image.
         *
         * @param destination
         * @throws IOException
         */
        private void toImage(File destination, boolean renderHeat) throws IOException {
            int[] rgba = new int[_w * 2 * _h * 4];

            //write ppm header
            try (BufferedWriter output = new BufferedWriter(new FileWriter(destination))) {
                //write ppm header
                output.write("P3\n" + (renderHeat ? _w * 2 : _w) + " " + _h + "\n" + 255 + "\n");

                for (int y = 0; y < _h; y++) {
                    for (int x = 0; x < _w; x++) {
                        int idxl = 0, idxr;
                        if (renderHeat) {
                            idxl = 4 * (x + y * _w * 2);
                            idxr = 4 * (x + y * _w * 2 + _w);
                        } else {
                            idxr = 4 * (x + y * _w);
                        }

                        double volume = _d.volume(x, y);

                        double shade = (1.0 - _d.getSrc(x, y)) * volume;
                        shade = min(max(shade, 0.0), 1.0);
                        rgba[idxr + 0] = (int) (shade * 255.0);
                        rgba[idxr + 1] = (int) (shade * 255.0);
                        rgba[idxr + 2] = (int) (shade * 255.0);
                        rgba[idxr + 3] = 0xFF;

                        if (_d._cell[x + y * _w] == CELL_EMPTY) {
                            rgba[idxr] = 0xFF;
                            rgba[idxr + 1] = rgba[idxr + 2] = 0;
                        }

                        if (renderHeat) {
                            double t = abs((float) _t.getSrc(x, y) - _tAmb) / 70.0;

                            t = min(max(t, 0.0), 1.0);

                            double r = 1.0 + volume * (min(t * 4.0, 1.0) - 1.0);
                            double g = 1.0 + volume * (min(t * 2.0, 1.0) - 1.0);
                            double b = 1.0 + volume * (max(min(t * 4.0 - 3.0, 1.0), 0.0) - 1.0);

                            rgba[idxl + 0] = (int) (r * 255.0);
                            rgba[idxl + 1] = (int) (g * 255.0);
                            rgba[idxl + 2] = (int) (b * 255.0);
                            rgba[idxl + 3] = 0xFF;
                        }
                    }
                }

                for (int i = 0; i < (renderHeat ? _w * 2 : _w) * _h; i++) {
                    output.write(rgba[i * 4] + " " + rgba[i * 4 + 1] + " " + rgba[i * 4 + 2] + " ");
                    if (i % _w == 0) {
                        output.write("\n");
                    }
                }
            }
        }
    }

    /**
     * Main class processing fluid particles
     */
    static class ParticleQuantities {

        /**
         * Maximum allowed number of particles per cell
         */
        static int _MaxPerCell = 12;

        /**
         * Minimum allowed number of particles per cell
         */
        static int _MinPerCell = 3;

        /**
         * Initial number of particles per cell
         */
        static int _AvgPerCell = 4;

        /**
         * Number of particles currently active
         */
        int _particleCount;

        /**
         * Maximum number of particles the simulation can handle
         */
        int _maxParticles;

        /**
         * The usual culprits
         */
        int _w;
        int _h;
        double _hx;

        /**
         * List of solid bodies to consider in the simulation
         */
        ArrayList<A_SolidBody> _bodies;

        /**
         * Filter weights (auxiliary array provided to fluid quantities)
         */
        double[] _weight;

        /**
         * Number of particles per cell
         */
        int[] _counts;

        /**
         * Particle positions
         */
        double[] _posX;
        double[] _posY;

        /**
         * Particle 'properties', that is, value for each fluid quantity
         * (velocities, density etc.)
         */
        ArrayList<double[]> _properties = new ArrayList<>();
        ArrayList<FluidQuantity> _quantities = new ArrayList<>();

        ParticleQuantities(
                int w, int h, double hx,
                ArrayList<A_SolidBody> bodies) {
            _w = w;
            _h = h;
            _hx = hx;
            _bodies = bodies;

            _maxParticles = _w * _h * _MaxPerCell;

            _posX = new double[_maxParticles];
            _posY = new double[_maxParticles];

            _weight = new double[(_w + 1) * (_h + 1)];
            _counts = new int[_w * _h];

            initParticles();
        }

        /**
         * Helper function returning true if a position is inside a solid body
         */
        private boolean pointInBody(double x, double y) {
            return _bodies.stream()
                    .anyMatch(_body -> _body.distance(x * _hx, y * _hx) < 0.0);
        }

        /**
         * Initializes particle positions on randomly jittered grid locations
         */
        private void initParticles() {
            int idx = 0;
            for (int y = 0; y < _h; y++) {
                for (int x = 0; x < _w; x++) {
                    for (int i = 0; i < _AvgPerCell; i++, idx++) {
                        _posX[idx] = x + (float) random();
                        _posY[idx] = y + (float) random();

                        /**
                         * Discard particles landing inside solid bodies
                         */
                        if (pointInBody(_posX[idx], _posY[idx])) {
                            idx--;
                        }
                    }
                }
            }

            _particleCount = idx;
        }

        /**
         * Counts the number of particles per cell
         */
        private void countParticles() {
            _counts = new int[_w * _h];
            for (int i = 0; i < _particleCount; i++) {
                int ix = (int) _posX[i];
                int iy = (int) _posY[i];

                if (ix >= 0 && iy >= 0 && ix < _w && iy < _h) {
                    _counts[ix + iy * _w]++;
                }
            }
        }

        /**
         * Decimates particles in crowded cells
         */
        private void pruneParticles() {
            for (int i = 0; i < _particleCount; i++) {
                int ix = (int) _posX[i];
                int iy = (int) _posY[i];
                int idx = ix + iy * _w;

                if (ix < 0 && iy < 0 && ix >= _w && iy >= _h) {
                    continue;
                }

                if (_counts[idx] > _MaxPerCell) {
                    int j = --_particleCount;
                    _posX[i] = _posX[j];
                    _posY[i] = _posY[j];
                    for (int t = 0; t < _quantities.size(); t++) {
                        _properties.get(t)[i] = _properties.get(t)[j];
                    }

                    _counts[idx]--;
                    i--;
                }
            }
        }

        /**
         * Adds new particles in cells with dangerously little particles
         */
        private void seedParticles() {
            for (int y = 0, idx = 0; y < _h; y++) {
                for (int x = 0; x < _w; x++, idx++) {
                    for (int i = 0; i < _MinPerCell - _counts[idx]; i++) {
                        if (_particleCount == _maxParticles) {
                            return;
                        }

                        int j = _particleCount;

                        _posX[j] = x + (float) random();
                        _posY[j] = y + (float) random();

                        /**
                         * Reject particle if it lands inside a solid body
                         */
                        if (pointInBody(_posX[idx], _posY[idx])) {
                            continue;
                        }

                        /**
                         * Get current grid values
                         */
                        for (int t = 0; t < _quantities.size(); t++) {
                            _properties.get(t)[j] = _quantities.get(t).lerpGrid(_posX[j], _posY[j]);
                        }

                        _particleCount++;
                    }
                }
            }
        }

        /**
         * Pushes particle back into the fluid if they land inside solid bodies
         */
        private Vector2D backProject(double x, double y) {
            double d = 1e30;
            int closestBody = -1;
            for (int i = 0; i < _bodies.size(); i++) {
                double id = _bodies.get(i).distance(x * _hx, y * _hx);

                if (id < d) {
                    d = id;
                    closestBody = i;
                }
            }

            if (d < -1.0) {
                x *= _hx;
                y *= _hx;
                Vector2D result = _bodies.get(closestBody).closestSurfacePoint(x, y);
                x = result.getX();
                y = result.getY();
                double nx = 0, ny = 0;
                result = _bodies.get(closestBody).distanceNormal(nx, ny, x, y);
                nx = result.getX();
                ny = result.getY();
                x -= nx * _hx;
                y -= ny * _hx;
                x /= _hx;
                y /= _hx;
            }

            return new Vector2D(x, y);
        }

        /**
         * The same Runge Kutta interpolation routine as before - only now
         * forward in time instead of backwards.
         */
        private Vector2D rungeKutta3(
                double x, double y,
                double timestep,
                FluidQuantity u, FluidQuantity v) {
            double firstU = u.lerpGrid(x, y) / _hx;
            double firstV = v.lerpGrid(x, y) / _hx;

            double midX = x + 0.5 * timestep * firstU;
            double midY = y + 0.5 * timestep * firstV;

            double midU = u.lerpGrid(midX, midY) / _hx;
            double midV = v.lerpGrid(midX, midY) / _hx;

            double lastX = x + 0.75 * timestep * midU;
            double lastY = y + 0.75 * timestep * midV;

            double lastU = u.lerpGrid(lastX, lastY);
            double lastV = v.lerpGrid(lastX, lastY);

            x += timestep * ((2.0 / 9.0) * firstU + (3.0 / 9.0) * midU + (4.0 / 9.0) * lastU);
            y += timestep * ((2.0 / 9.0) * firstV + (3.0 / 9.0) * midV + (4.0 / 9.0) * lastV);

            return new Vector2D(x, y);
        }

        /**
         * Adds a new quantity to be carried by the particles
         */
        private void addQuantity(FluidQuantity q) {
            double[] property = new double[_maxParticles];

            _quantities.add(q);
            _properties.add(property);
        }

        /**
         * Interpolates the change in quantity back onto the particles. Mixes in
         * a little bit of the pure Particle-in-cell update using the parameter
         * alpha.
         */
        private void gridToParticles(double alpha) {
            for (int t = 0; t < _quantities.size(); t++) {
                for (int i = 0; i < _particleCount; i++) {
                    _properties.get(t)[i] *= 1.0 - alpha;
                    _properties.get(t)[i] += _quantities.get(t).lerpGrid(_posX[i], _posY[i]);
                }
            }
        }

        /**
         * Interpolates particle quantities onto the grid, extrapolates them and
         * spawns/prunes particles where necessary.
         */
        private void particlesToGrid() {
            for (int t = 0; t < _quantities.size(); t++) {
                _quantities.get(t).fromParticles(_weight, _particleCount, _posX, _posY, _properties.get(t));
                _quantities.get(t).extrapolate();
            }

            countParticles();
            pruneParticles();
            seedParticles();

            log.log(INFO, format("Particle count: %d", _particleCount));
        }

        /**
         * Advects particle in velocity field and clamps resulting positions to
         * the fluid domain
         */
        private void advect(double timestep, FluidQuantity u, FluidQuantity v) {
            for (int i = 0; i < _particleCount; i++) {
                Vector2D result = rungeKutta3(_posX[i], _posY[i], timestep, u, v);
                _posX[i] = result.getX();
                _posY[i] = result.getY();
                result = backProject(_posX[i], _posY[i]);
                _posX[i] = result.getX();
                _posY[i] = result.getY();
                _posX[i] = max(min(_posX[i], _w - 0.001), 0.0);
                _posY[i] = max(min(_posY[i], _h - 0.001), 0.0);
            }
        }
    }

    /**
     * Enum to differentiate fluid and solid cells
     */
    public enum CellType {

        CELL_FLUID,
        CELL_SOLID,
        CELL_EMPTY
    }
}
