#pragma once

//          Copyright Erik Lundin 2016.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

// Version: 1.0.0

#include <numbers>


namespace gausskruger {


struct Ellipsoid
{
  double flattening;
  double equatorialRadius;

  constexpr double e2() const
  {
    return flattening * (2 - flattening); // e2: first eccentricity squared
  }
  constexpr double n() const
  {
    return flattening / (2 - flattening); // n: 3rd flattening
  }
  constexpr double rectifyingRadius() const
  {
    return equatorialRadius / (1 + n()) * (1 + 0.25 * std::pow(n(), 2) + 0.015625 * std::pow(n(), 4));
  }
};

constexpr auto GRS80 = Ellipsoid{ .flattening = 1.0 / 298.257222101, .equatorialRadius = 6378137.0 };

struct Projection
{
  Ellipsoid ellipsoid;
  double centralMeridian;
  double scale;
  double falseNorthing;
  double falseEasting;
};

constexpr auto SWEREF99TM = Projection{ .ellipsoid = GRS80, .centralMeridian = 15.0, .scale = 0.9996, .falseNorthing = 0.0, .falseEasting = 500000.0 };

struct GridCoord
{
  double northing;
  double easting;
};

struct GeodeticCoord
{
  double latitude;
  double longitude;
};

template <Projection Projection>
constexpr GridCoord geodeticToGrid(GeodeticCoord coord)
{
    const double e2 = Projection.ellipsoid.e2();
    const double n = Projection.ellipsoid.n();
    const double rectifyingRadius = Projection.ellipsoid.rectifyingRadius();

    double A = e2;
    double B = (5 * std::pow(e2, 2) - std::pow(e2, 3)) / 6.0;
    double C = (104 * std::pow(e2, 3) - 45 * std::pow(e2, 4)) / 120.0;
    double D = (1237 * std::pow(e2, 4)) / 1260.0;

    // Latitude and longitude are expected to be given in degrees
    // phi and lambda: latitude and longitude in radians
    double phi = coord.latitude * M_PI / 180;
    double lambda = coord.longitude * M_PI / 180;
    double lambda0 = Projection.centralMeridian * M_PI / 180;

    // deltaLambda: longitude relative to the central meridian
    double deltaLambda = lambda - lambda0;

    // phiStar: conformal latitude
    double phiStar =
            phi - sin(phi) * cos(phi) *
            (A + B*std::pow(sin(phi), 2) + C*std::pow(sin(phi), 4) + D*std::pow(sin(phi), 6));

    double xiPrim = atan(tan(phiStar) / cos(deltaLambda));
    double etaPrim = atanh(cos(phiStar) * sin(deltaLambda));

    double beta1 = 1/2.0 * n - 2/3.0 * std::pow(n, 2) + 5/16.0 * std::pow(n, 3)     + 41/180.0 * std::pow(n, 4);
    double beta2 =           13/48.0 * std::pow(n, 2)  - 3/5.0 * std::pow(n, 3)   + 557/1440.0 * std::pow(n, 4);
    double beta3 =                               61/240.0 * std::pow(n, 3)    - 103/140.0 * std::pow(n, 4);
    double beta4 =                                                    49561/161280.0 * std::pow(n, 4);

    auto northing = Projection.falseNorthing
            + Projection.scale * rectifyingRadius * (xiPrim
                                            + beta1 * sin(2*xiPrim) * cosh(2*etaPrim)
                                            + beta2 * sin(4*xiPrim) * cosh(4*etaPrim)
                                            + beta3 * sin(6*xiPrim) * cosh(6*etaPrim)
                                            + beta4 * sin(8*xiPrim) * cosh(8*etaPrim));
    auto easting = Projection.falseEasting
            + Projection.scale * rectifyingRadius * (etaPrim
                                            + beta1 * cos(2*xiPrim) * sinh(2*etaPrim)
                                            + beta2 * cos(4*xiPrim) * sinh(4*etaPrim)
                                            + beta3 * cos(6*xiPrim) * sinh(6*etaPrim)
                                            + beta4 * cos(8*xiPrim) * sinh(8*etaPrim));
    return { northing, easting };
}

template <Projection Projection>
constexpr GeodeticCoord gridToGeodetic(GridCoord coord)
{
    const double e2 = Projection.ellipsoid.e2();
    const double n = Projection.ellipsoid.n();
    const double rectifyingRadius = Projection.ellipsoid.rectifyingRadius();

    double xi = (coord.northing - Projection.falseNorthing) / (Projection.scale * rectifyingRadius);
    double eta = (coord.easting - Projection.falseEasting) / (Projection.scale * rectifyingRadius);

    double delta1 = 1/2.0 * n - 2/3.0 * std::pow(n, 2) + 37/96.0 * std::pow(n, 3)     - 1/360.0 * std::pow(n, 4);
    double delta2 =            1/48.0 * std::pow(n, 2)  + 1/15.0 * std::pow(n, 3)  - 437/1440.0 * std::pow(n, 4);
    double delta3 =                                17/480.0 * std::pow(n, 3)    - 37/840.0 * std::pow(n, 4);
    double delta4 =                                                     4397/161280.0 * std::pow(n, 4);

    double xiPrim = xi
            - delta1 * sin(2*xi) * cosh(2*eta)
            - delta2 * sin(4*xi) * cosh(4*eta)
            - delta3 * sin(6*xi) * cosh(6*eta)
            - delta4 * sin(8*xi) * cosh(8*eta);
    double etaPrim = eta
            - delta1 * cos(2*xi) * sinh(2*eta)
            - delta2 * cos(4*xi) * sinh(4*eta)
            - delta3 * cos(6*xi) * sinh(6*eta)
            - delta4 * cos(8*xi) * sinh(8*eta);

    double phiStar = asin(sin(xiPrim) / cosh(etaPrim)); // Conformal latitude
    double deltaLambda = atan(sinh(etaPrim) / cos(xiPrim));

    double AStar =  e2     + std::pow(e2, 2)       + std::pow(e2, 3)        + std::pow(e2, 4);
    double BStar =      (7 * std::pow(e2, 2)  + 17 * std::pow(e2, 3)   + 30 * std::pow(e2, 4)) / -6;
    double CStar =                       (224 * std::pow(e2, 3)  + 889 * std::pow(e2, 4)) / 120;
    double DStar =                                          (4279 * std::pow(e2, 4)) / -1260;

    double phi = phiStar
            + sin(phiStar) * cos(phiStar) * (  AStar
                                             + BStar * std::pow(sin(phiStar), 2)
                                             + CStar * std::pow(sin(phiStar), 4)
                                             + DStar * std::pow(sin(phiStar), 6));

    // phi: latitude in radians, lambda: longitude in radians
    // Return latitude and longitude as degrees
    return {
      .latitude = phi * 180.0 / std::numbers::pi,
      .longitude = Projection.centralMeridian + deltaLambda * 180.0 / std::numbers::pi
    };
}

} // namespace gausskruger
