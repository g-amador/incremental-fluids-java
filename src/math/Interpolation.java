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
package math;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.util.logging.Level.OFF;
import java.util.logging.Logger;
import static java.util.logging.Logger.getLogger;

/**
 * Class that implements numerical interpolation algorithms.
 *
 * @author G. Amador {@literal &} A. Gomes
 */
public class Interpolation {

    static final Logger log = getLogger("Interpolation");

    static {
        log.setLevel(OFF);
    }

    /**
     * Single precision linear interpolate between a and b for x ranging from 0
     * to 1.
     *
     * @param a
     * @param b
     * @param x
     * @return single precision linear interpolation between a and b for x
     * ranging from 0 to 1.
     */
    public static float lerp(
            float a, float b,
            float x) {
        // Imprecise method which does not guarantee a  = b when x = 1,
        // due to floating-point arithmetic error.
        // return a + x * (b - a);

        // Precise method which guarantees a = b when x = 1.
        return (1 - x) * a + x * b;
    }

    /**
     * Double precision linear interpolate between a and b for x ranging from 0
     * to 1.
     *
     * @param a
     * @param b
     * @param x
     * @return double precision linear interpolation between a and b for x
     * ranging from 0 to 1.
     */
    public static double lerp(
            double a, double b,
            double x) {
        // Imprecise method which does not guarantee a  = b when x = 1,
        // due to floating-point arithmetic error.
        // return a + x * (b - a);

        // Precise method which guarantees a = b when x = 1.
        return (1 - x) * a + x * b;
    }

    /**
     * Single precision cubic interpolate using samples a through d for x
     * ranging from 0 to 1. A Catmull-Rom spline is used. Over- and undershoots
     * are clamped to prevent blow-up.
     *
     * @param a
     * @param b
     * @param c
     * @param d
     * @param x
     * @return single precision cubic interpolation using samples a through d
     * for x ranging from 0 to 1. A Catmull-Rom spline is used. Over- and
     * undershoots are clamped to prevent blow-up.
     */
    public static float cerp(
            float a, float b,
            float c, float d,
            float x) {
        float xsq = x * x;
        float xcu = xsq * x;

        float minV = min(a, min(b, min(c, d)));
        float maxV = max(a, max(b, max(c, d)));

        float t
                = a * (0.0f - 0.5f * x + 1.0f * xsq - 0.5f * xcu)
                + b * (1.0f + 0.0f * x - 2.5f * xsq + 1.5f * xcu)
                + c * (0.0f + 0.5f * x + 2.0f * xsq - 1.5f * xcu)
                + d * (0.0f + 0.0f * x - 0.5f * xsq + 0.5f * xcu);

        return min(max(t, minV), maxV);
    }

    /**
     * Double precision cubic interpolate using samples a through d for x
     * ranging from 0 to 1. A Catmull-Rom spline is used. Over- and undershoots
     * are clamped to prevent blow-up.
     *
     * @param a
     * @param b
     * @param c
     * @param d
     * @param x
     * @return double precision cubic interpolation using samples a through d
     * for x ranging from 0 to 1. A Catmull-Rom spline is used. Over- and
     * undershoots are clamped to prevent blow-up.
     */
    public static double cerp(
            double a, double b,
            double c, double d,
            double x) {
        double xsq = x * x;
        double xcu = xsq * x;

        double minV = min(a, min(b, min(c, d)));
        double maxV = max(a, max(b, max(c, d)));

        double t
                = a * (0.0 - 0.5 * x + 1.0 * xsq - 0.5 * xcu)
                + b * (1.0 + 0.0 * x - 2.5 * xsq + 1.5 * xcu)
                + c * (0.0 + 0.5 * x + 2.0 * xsq - 1.5 * xcu)
                + d * (0.0 + 0.0 * x - 0.5 * xsq + 0.5 * xcu);

        return min(max(t, minV), maxV);
    }
}
