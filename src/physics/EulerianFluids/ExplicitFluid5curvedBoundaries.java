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
import static java.lang.Double.isNaN;
import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
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
import static math.Interpolation.cerp;
import static math.Interpolation.lerp;
import org.apache.commons.math3.geometry.euclidean.twod.Vector2D;
import static physics.EulerianFluids.ExplicitFluid5curvedBoundaries.CellType.CELL_FLUID;
import static physics.EulerianFluids.ExplicitFluid5curvedBoundaries.CellType.CELL_SOLID;
import physics.solidBodies.A_SolidBody;
import physics.solidBodies.SolidSphere;

/**
 * @author Gonçalo Amador {@literal &} Abel Gomes
 */
class ExplicitFluid5curvedBoundaries {

    static final Logger log = getLogger("ExplicitFluid5curvedBoundaries");

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

        ArrayList<A_SolidBody> bodies = new ArrayList<>();
        //bodies.add(new SolidBox(0.5, 0.6, 0.7, 0.1, PI * 0.25, 0.0, 0.0, 0.0));
        bodies.add(new SolidSphere(0.5, 0.6, 0.1, PI * 0.25, 0.0, 0.0, 0.0));

        FluidSolver solver = new FluidSolver(sizeX, sizeY, density, bodies);

        double time = 0.0;
        int iterations = 0;

        while (time < 8.0) {
            for (int i = 0; i < 4; i++) {
                solver.addInflow(0.45, 0.2, 0.15, 0.03, 1.0, 0.0, 3.0);
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

            for (A_SolidBody body : bodies) {
                body.update(timestep);
            }
        }
    }

    /**
     * Length of vector (x, y)
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
            /**
             * All outside
             */
            case 0x0:
                return 0.0;
            /**
             * One inside
             */
            case 0x1:
                return triangleOccupancy(d21, d11, d12);
            case 0x2:
                return triangleOccupancy(d11, d12, d22);
            case 0x4:
                return triangleOccupancy(d12, d22, d21);
            case 0x8:
                return triangleOccupancy(d22, d21, d11);
            /**
             * One outside
             */
            case 0xE:
                return 1.0 - triangleOccupancy(-d21, -d11, -d12);
            case 0xD:
                return 1.0 - triangleOccupancy(-d11, -d12, -d22);
            case 0xB:
                return 1.0 - triangleOccupancy(-d12, -d22, -d21);
            case 0x7:
                return 1.0 - triangleOccupancy(-d22, -d21, -d11);
            /**
             * Two adjacent inside
             */
            case 0x3:
                return trapezoidOccupancy(d21, d22, d11, d12);
            case 0x6:
                return trapezoidOccupancy(d11, d21, d12, d22);
            case 0x9:
                return trapezoidOccupancy(d12, d22, d11, d21);
            case 0xC:
                return trapezoidOccupancy(d11, d12, d21, d22);
            /**
             * Two opposed inside
             */
            case 0x5:
                return triangleOccupancy(d11, d12, d22)
                        + triangleOccupancy(d22, d21, d11);
            case 0xA:
                return triangleOccupancy(d21, d11, d12)
                        + triangleOccupancy(d12, d22, d21);
            /**
             * All inside
             */
            case 0xF:
                return 1.0;
        }

        return 0.0;
    }

    static class FluidQuantity {

        /**
         * Memory buffers for fluid quantity
         */
        double[] _src;
        double[] _dst;
        double[] _temp;

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
            _dst = new double[_w * _h];

            /* Make distance grid one larger in each dimension */
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

        private double volume(int x, int y) {
            return _volume[x + y * _w];
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
         * Cubic interpolate on grid at coordinates (x, y). Coordinates will be
         * clamped to lie in simulation domain
         */
        private double cerpGrid(double x, double y) {
            x = min(max(x - _ox, 0.0), _w - 1.001);
            y = min(max(y - _oy, 0.0), _h - 1.001);
            int ix = (int) x;
            int iy = (int) y;

            x -= ix;
            y -= iy;

            int x0 = max(ix - 1, 0), x1 = ix, x2 = ix + 1, x3 = min(ix + 2, _w - 1);
            int y0 = max(iy - 1, 0), y1 = iy, y2 = iy + 1, y3 = min(iy + 2, _h - 1);

            double q0 = cerp(getSrc(x0, y0), getSrc(x1, y0), getSrc(x2, y0), getSrc(x3, y0), x);
            double q1 = cerp(getSrc(x0, y1), getSrc(x1, y1), getSrc(x2, y1), getSrc(x3, y1), x);
            double q2 = cerp(getSrc(x0, y2), getSrc(x1, y2), getSrc(x2, y2), getSrc(x3, y2), x);
            double q3 = cerp(getSrc(x0, y3), getSrc(x1, y3), getSrc(x2, y3), getSrc(x3, y3), x);

            return cerp(q0, q1, q2, q3, y);
        }

        /**
         * If the point (x, y) is inside a solid, project it back out to the
         * closest point on the surface of the solid.
         */
        private Vector2D backProject(double x, double y, ArrayList<A_SolidBody> bodies) {
            int rx = min(max((int) (x - _ox), 0), _w - 1);
            int ry = min(max((int) (y - _oy), 0), _h - 1);

            if (_cell[rx + ry * _w] != CELL_FLUID) {
                x = (x - _ox) * _hx;
                y = (y - _oy) * _hx;

                Vector2D result = bodies.get(_body[rx + ry * _w]).closestSurfacePoint(x, y);
                x = result.getX() / _hx + _ox;
                y = result.getY() / _hx + _oy;
            }

            return new Vector2D(x, y);
        }

        /**
         * Advect grid in velocity field u, v with given timestep
         */
        private void advect(double timestep, FluidQuantity u, FluidQuantity v,
                ArrayList<A_SolidBody> bodies) {

            for (int iy = 0, idx = 0; iy < _h; iy++) {
                for (int ix = 0; ix < _w; ix++, idx++) {
                    if (_cell[idx] == CELL_FLUID) {
                        double x = ix + _ox;
                        double y = iy + _oy;

                        /**
                         * First component: Integrate in time
                         */
                        Vector2D result = RungeKutta3(x, y, timestep, u, v);
                        x = result.getX();
                        y = result.getY();

                        /**
                         * If integrating back in time leaves us inside a solid
                         * boundary (due to numerical error), make sure we
                         * interpolate from a point inside the fluid.
                         */
                        result = backProject(x, y, bodies);
                        x = result.getX();
                        y = result.getY();

                        /**
                         * Second component: Interpolate from grid
                         */
                        _dst[idx] = cerpGrid(x, y);
                    }
                }
            }
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

                    Vector2D result = bodies.get(_body[idx]).distanceNormal(_normalX[idx], _normalY[idx], x, y);
                    _normalX[idx] = result.getX();
                    _normalY[idx] = result.getY();

                    /**
                     * Solid cells are now defined as cells with zero fluid
                     * volume
                     */
                    if (_volume[idx] == 0.0) {
                        _cell[idx] = CELL_SOLID;
                    } else {
                        _cell[idx] = CELL_FLUID;
                    }
                }
            }
        }

        /**
         * Prepare auxiliary array for extrapolation. The purpose of
         * extrapolation is to extrapolate fluid quantities into solids, where
         * these quantities would normally be undefined. However, we need these
         * values for stable interpolation and boundary conditions.
         *
         * The way these are extrapolated here is by essentially solving a PDE,
         * such that the gradient of the fluid quantity is 0 along the gradient
         * of the distance field. This is essentially a more robust formulation
         * of "Set quantity inside solid to the value at the closest point on
         * the solid-fluid boundary"
         *
         * This PDE has a particular form which makes it very easy to solve
         * exactly using an upwinding scheme. What this means is that we can
         * solve it from outside-to-inside, with information flowing along the
         * normal from the boundary.
         *
         * Specifically, we can solve for the value inside a fluid cell using
         * extrapolateNormal if the two adjacent grid cells in "upstream"
         * direction (where the normal points to) are either fluid cells or have
         * been solved for already.
         *
         * The mask array keeps track of which cells wait for which values. If
         * an entry is 0, it means both neighbours are available and the cell is
         * ready for the PDE solve. If it is 1, the cell waits for the neighbour
         * in x direction, 2 for y-direction and 3 for both.
         */
        private void fillSolidMask() {
            for (int y = 1; y < _h - 1; y++) {
                for (int x = 1; x < _w - 1; x++) {
                    int idx = x + y * _w;

                    if (_cell[idx] == CELL_FLUID) {
                        continue;
                    }

                    double nx = _normalX[idx];
                    double ny = _normalY[idx];

                    _mask[idx] = 0;
                    if (nx != 0.0 && _cell[idx + (int) signum(nx)] != CELL_FLUID) {
                        _mask[idx] |= 1; /* Neighbour in normal x direction is blocked */

                    }
                    if (ny != 0.0 && _cell[idx + (int) signum(ny) * _w] != CELL_FLUID) {
                        _mask[idx] |= 2; /* Neighbour in normal y direction is blocked */

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
         * Given that a neighbour in upstream direction specified by mask (1=x,
         * 2=y) now has been solved for, update the mask appropriately and, if
         * this cell can now be computed, add it to the queue of ready cells
         */
        private Stack<Integer> freeNeighbour(
                int idx,
                Stack<Integer> border,
                int mask) {
            _mask[idx] &= ~mask;
            if (_cell[idx] != CELL_FLUID && _mask[idx] == 0) {
                border.push(idx);
            }
            return border;
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

                /**
                 * Solve for value in cell
                 */
                _src[idx] = extrapolateNormal(idx);

                /**
                 * Notify adjacent cells that this cell has been computed and
                 * can be used as an upstream value
                 */
                if (_normalX[idx - 1] > 0.0) {
                    border = freeNeighbour(idx - 1, border, 1);
                }
                if (_normalX[idx + 1] < 0.0) {
                    border = freeNeighbour(idx + 1, border, 1);
                }
                if (_normalY[idx - _w] > 0.0) {
                    border = freeNeighbour(idx - _w, border, 2);
                }
                if (_normalY[idx + _w] < 0.0) {
                    border = freeNeighbour(idx + _w, border, 2);
                }
            }
        }

        /**
         * Third order Runge-Kutta for velocity integration in time
         */
        private Vector2D RungeKutta3(double x, double y, double timestep, FluidQuantity u, FluidQuantity v) {
            double firstU = u.lerpGrid(x, y) / _hx;
            double firstV = v.lerpGrid(x, y) / _hx;

            double midX = x - 0.5 * timestep * firstU;
            double midY = y - 0.5 * timestep * firstV;

            double midU = u.lerpGrid(midX, midY) / _hx;
            double midV = v.lerpGrid(midX, midY) / _hx;

            double lastX = x - 0.75 * timestep * midU;
            double lastY = y - 0.75 * timestep * midV;

            double lastU = u.lerpGrid(lastX, lastY);
            double lastV = v.lerpGrid(lastX, lastY);

            x -= timestep * ((2.0 / 9.0) * firstU + (3.0 / 9.0) * midU + (4.0 / 9.0) * lastU);
            y -= timestep * ((2.0 / 9.0) * firstV + (3.0 / 9.0) * midV + (4.0 / 9.0) * lastV);

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
         * List of solid bodies to consider in the simulation
         */
        ArrayList<A_SolidBody> _bodies;

        FluidSolver(int w, int h, double density, ArrayList<A_SolidBody> bodies) {
            _w = w;
            _h = h;
            _density = density;
            _bodies = bodies;

            _hx = 1.0 / min(w, h);

            _d = new FluidQuantity(_w, _h, 0.5, 0.5, _hx);
            _u = new FluidQuantity(_w + 1, _h, 0.0, 0.5, _hx);
            _v = new FluidQuantity(_w, _h + 1, 0.5, 0.0, _hx);

            _r = new double[_w * _h];
            _p = new double[_w * _h];
            _z = new double[_w * _h];
            _s = new double[_w * _h];

            _aDiag = new double[_w * _h];
            _aPlusX = new double[_w * _h];
            _aPlusY = new double[_w * _h];
            _precon = new double[_w * _h];
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
         * Entries of the pressure matrix are modified accordingly
         */
        private void buildPressureMatrix(double timestep) {
            double scale = timestep / (_density * _hx * _hx);
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
                        double factor = scale * _u.volume(x + 1, y);
                        _aDiag[idx] += factor;
                        _aDiag[idx + 1] += factor;
                        _aPlusX[idx] = -factor;
                    }
                    if (y < _h - 1 && cell[idx + _w] == CELL_FLUID) {
                        double factor = scale * _v.volume(x, y + 1);
                        _aDiag[idx] += factor;
                        _aDiag[idx + _w] += factor;
                        _aPlusY[idx] = -factor;
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
            double result = 0.0;
            for (int i = 0; i < _w * _h; i++) {
                result += a[i] * b[i];
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
            for (int i = 0; i < _w * _h; i++) {
                dst[i] = a[i] + b[i] * s;
            }
        }

        /**
         * Returns maximum absolute value in vector 'a'
         */
        private double infinityNorm(double[] a) {
            double maxA = 0.0;
            for (int i = 0; i < _w * _h; i++) {
                maxA = max(maxA, abs((float) a[i]));
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
         * Applies the computed pressure to the velocity field
         */
        private void applyPressure(double timestep) {
            double scale = timestep / (_density * _hx);
            CellType[] cell = _d._cell;

            for (int y = 0, idx = 0; y < _h; y++) {
                for (int x = 0; x < _w; x++, idx++) {
                    if (cell[idx] != CELL_FLUID) {
                        continue;
                    }

                    _u.setSrc(x, y, _u.getSrc(x, y) - scale * _p[idx]);
                    _u.setSrc(x + 1, y, _u.getSrc(x + 1, y) + scale * _p[idx]);
                    _v.setSrc(x, y, _v.getSrc(x, y) - scale * _p[idx]);
                    _v.setSrc(x, y + 1, _v.getSrc(x, y + 1) + scale * _p[idx]);
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
            _u.fillSolidFields(_bodies);
            _v.fillSolidFields(_bodies);

            setBoundaryCondition();

            buildRhs();
            buildPressureMatrix(timestep);
            buildPreconditioner();
            project(2_000);
            applyPressure(timestep);

            _d.extrapolate();
            _u.extrapolate();
            _v.extrapolate();

            setBoundaryCondition();

            _d.advect(timestep, _u, _v, _bodies);
            _u.advect(timestep, _u, _v, _bodies);
            _v.advect(timestep, _u, _v, _bodies);

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

                for (int i = 0; i < _w * _h; i++) {
                    /* Use fluid volume for nice anti aliasing */
                    int shade = (int) ((1.0 - _d._src[i]) * _d.volume(i % _w, i / _w) * 255.0);
                    shade = max(min(shade, 255), 0);

                    output.write(shade + " " + shade + " " + shade + " ");
                    if (i % _w == 0) {
                        output.write("\n");
                    }
                }
            }
        }
    }

    /**
     * Enum to differentiate fluid and solid cells
     */
    public enum CellType {

        CELL_FLUID,
        CELL_SOLID
    }
}
