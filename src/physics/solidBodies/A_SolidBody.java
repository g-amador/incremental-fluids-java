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
package physics.solidBodies;

import static java.lang.Math.cos;
import static java.lang.Math.sin;
import static java.lang.Math.sqrt;
import org.apache.commons.math3.geometry.euclidean.twod.Vector2D;

/**
 * The base abstract class representing solid bodies in the simulation. It holds
 * information about position, scale and rotation of the solid as well as
 * lateral and angular velocity.
 *
 * It does not represent any shape; this is handled by the subclasses. To expose
 * the shape to the simulation, methods for evaluating the signed distance to
 * the solid, the gradient of the distance function and the closest point on the
 * surface of the solid are exposed.
 *
 * @author G. Amador {@literal &} A. Gomes
 */
public abstract class A_SolidBody {

    /**
     * Length of ArrayList (x, y)
     *
     * @param x
     * @param y
     * @return
     */
    static protected double length(double x, double y) {
        return sqrt(x * x + y * y);
    }

    /**
     * Rotates point (x, y) by angle phi
     *
     * @param x
     * @param y
     * @param phi
     * @return
     */
    static protected Vector2D rotate(double x, double y, double phi) {
        double tmpX = x, tmpY = y;
        x = cos(phi) * tmpX + sin(phi) * tmpY;
        y = -sin(phi) * tmpX + cos(phi) * tmpY;

        return new Vector2D(x, y);
    }

    protected double _posX; /* Position */

    protected double _posY;

    protected double _scaleX; /* Scale */

    protected double _scaleY;

    protected double _theta; /* Rotation */

    protected double _velX; /* Lateral velocity */

    protected double _velY;

    protected double _velTheta; /* Angular velocity */


    /**
     * Constructor.
     *
     * @param posX
     * @param posY
     * @param scaleX
     * @param scaleY
     * @param theta
     * @param velX
     * @param velY
     * @param velTheta
     */
    A_SolidBody(double posX, double posY, double scaleX, double scaleY,
            double theta, double velX, double velY, double velTheta) {
        _posX = posX;
        _posY = posY;
        _scaleX = scaleX;
        _scaleY = scaleY;
        _theta = theta;
        _velX = velX;
        _velY = velY;
        _velTheta = velTheta;
    }

    /**
     * Transforms point (x, y) form the global to the local coordinate system
     *
     * @param x
     * @param y
     * @return
     */
    protected Vector2D globalToLocal(double x, double y) {
        x -= _posX;
        y -= _posY;
        Vector2D result = rotate(x, y, -_theta);
        x = result.getX() / _scaleX;
        y = result.getY() / _scaleY;

        return new Vector2D(x, y);
    }

    /**
     * Transforms point (x, y) form the local to the global coordinate system
     *
     * @param x
     * @param y
     * @return
     */
    protected Vector2D localToGlobal(double x, double y) {
        x *= _scaleX;
        y *= _scaleY;
        Vector2D result = rotate(x, y, _theta);
        x = result.getX() + _posX;
        y = result.getY() + _posY;

        return new Vector2D(x, y);
    }

    /**
     * Returns the signed distance from (x, y) to the nearest point on surface
     * of the solid. The distance is negative if (x, y) is inside the solid
     *
     * @param x
     * @param y
     * @return
     */
    public abstract double distance(double x, double y);
//    public double distance(double x, double y) {
//        return 0;
//    }

    /**
     * Changes (x, y) to lie on the closest point on the surface of the solid
     *
     * @param x
     * @param y
     * @return
     */
    public abstract Vector2D closestSurfacePoint(double x, double y);
//    public Vector2D closestSurfacePoint(double x, double y) {
//        return new Vector2D(x, y);
//    }

    /**
     * Returns the gradient of the distance function at (x, y) in (nx, ny)
     *
     * @param nx
     * @param ny
     * @param x
     * @param y
     * @return
     */
    public abstract Vector2D distanceNormal(double nx, double ny, double x, double y);
//    public Vector2D distanceNormal(double nx, double ny, double x, double y) {
//        return new Vector2D(nx, ny);
//    }

    /**
     * Evaluates velocities of the solid at a given point
     *
     * @param x
     * @param y
     * @return
     */
    public double velocityX(double x, double y) {
        return (_posY - y) * _velTheta + _velX;
    }

    /**
     * Evaluates velocities of the solid at a given point
     *
     * @param x
     * @param y
     * @return
     */
    public double velocityY(double x, double y) {
        return (x - _posX) * _velTheta + _velY;
    }

    /**
     * Simple Euler integration - enough for solid bodies, since they are not
     * influenced by the simulation and velocities are typically static
     *
     * @param timestep
     */
    public void update(double timestep) {
        _posX += _velX * timestep;
        _posY += _velY * timestep;
        _theta += _velTheta * timestep;
    }
}
