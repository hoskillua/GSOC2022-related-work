/**
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 **/

 /**
  * @file geometry/meshes/curvature-measures-icnc-3d.cpp
  * @ingroup Examples
  * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
  * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
  *
  * @date 2021/10/25
  *
  * An example file named curvature-measures-icnc-3d.
  *
  * This file is part of the DGtal library.
  */

  /**
     Computation of curvatures on a torus mesh, using interpolated
     corrected curvature measures (based on the theory of corrected
     normal currents).

  \verbatim
  ./examples/geometry/meshes/curvature-measures-icnc-3d torus 20 20 0.5
  \endverbatim
  outputs
  \verbatim
  Expected mean curvatures: min=0.25 max=0.625
  Computed mean curvatures: min=0.264763 max=0.622318
  Expected Gaussian curvatures: min=-0.5 max=0.25
  Computed Gaussian curvatures: min=-0.470473 max=0.244636
  \endverbatim

  It also produces several OBJ files to display curvature estimation
  results, `example-cnc-H.obj` and `example-cnc-G.obj` as well as the
  associated MTL file.

  <table>
  <tr><td>
  \image html torus-cnc-H-True-r0.jpg "Interpolated corrected mean curvature measure, r=0" width=90%
  </td><td>
  \image html torus-cnc-G-True-r0.jpg "Interpolated corrected Gaussian curvature measure, r=0" width=90%
  </td></tr>
  <tr><td>
  \image html torus-cnc-H-True-r0_5.jpg "Interpolated corrected mean curvature measure, r=0.5" width=90%
  </td><td>
  \image html torus-cnc-G-True-r0_5.jpg "Interpolated corrected Gaussian curvature measure, r=0.5" width=90%
  </td></tr>
  </table>

  @note Interpolated corected curvature measures can provide consistent results even on on bad sampling of smooth surfaces. A well known example is the Schwarz lantern. You may try the following line:

  \verbatim
  ./examples/geometry/meshes/curvature-measures-icnc-3d lantern 20 20 0.5
  \endverbatim

  outputs

  \verbatim
  Expected mean curvatures: min=0.25 max=0.25
  Computed mean curvatures: min=0.25 max=0.25
  Expected Gaussian curvatures: min=0 max=0
  Computed Gaussian curvatures: min=0 max=0
  \endverbatim

  @see \ref moduleCurvatureMeasures


  \example geometry/meshes/curvature-measures-icnc-3d.cpp
  */

#include <iostream>
#include <algorithm>
#include "DGtal/base/Common.h"
#include "DGtal/shapes/SurfaceMesh.h"
#include "DGtal/shapes/SurfaceMeshHelper.h"
  //! [curvature-measures-Includes]
#include "DGtal/geometry/meshes/CorrectedNormalCurrentComputer.h"
//! [curvature-measures-Includes]
#include "DGtal/io/writers/SurfaceMeshWriter.h"
#include "DGtal/io/readers/SurfaceMeshReader.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/io/colormaps/QuantifiedColorMap.h"
#include <cstdlib>
#include <chrono>

DGtal::GradientColorMap< double >
makeColorMap(double min_value, double max_value)
{
    DGtal::GradientColorMap< double > gradcmap(min_value, max_value);
    gradcmap.addColor(DGtal::Color(0, 0, 255));
    gradcmap.addColor(DGtal::Color(0, 255, 255));
    gradcmap.addColor(DGtal::Color(255, 255, 255));
    gradcmap.addColor(DGtal::Color(255, 255, 0));
    gradcmap.addColor(DGtal::Color(255, 0, 0));
    return gradcmap;
}

void usage(int argc, char* argv[])
{
    std::cout << "Usage: " << std::endl
        << "\t" << argv[0] << " <shape> <m> <n> <R>" << std::endl
        << std::endl
        << "Computation of mean and Gaussian curvatures on a shape, " << std::endl
        << "using interpolated corrected curvature measures (based " << std::endl
        << "on the theory of corrected normal currents)." << std::endl
        << "- builds a <shape> in {torus,lantern,sphere}, with     " << std::endl
        << "  <m> latitude points and <n> longitude points." << std::endl
        << "- <R> is the radius of the measuring balls." << std::endl
        << "It produces several OBJ files to display mean and" << std::endl
        << "Gaussian curvature estimation results: `example-cnc-H.obj`" << std::endl
        << "and `example-cnc-G.obj` as well as the associated MTL file." << std::endl;
}

int main(int argc, char* argv[])
{
    std::vector<std::string> O = { "torus","sphere","lantern" };
    std::vector<double> M = { 10,20,50, 75,100 };
    std::vector<double> N = { 10,20,50, 75,100 };
    std::vector<double> R = { 0,0.1,0.5 };

    //! [curvature-measures-Typedefs]
    using namespace DGtal;
    using namespace DGtal::Z3i;
    typedef SurfaceMesh< RealPoint, RealVector >                    SM;
    typedef CorrectedNormalCurrentComputer< RealPoint, RealVector > CNC;
    typedef SurfaceMeshHelper< RealPoint, RealVector >              SMH;

    long long quad_time_0 = 0;
    long long tri_time_0 = 0;
    long long quad_time_1 = 0;
    long long tri_time_1 = 0;
    long long quad_time_2 = 0;
    long long tri_time_2 = 0;

    for (auto input : O)
        for (auto r : R)
        {


           double quad_err_mean = 0;
           double tri_err_mean = 0;
           double quad_err_gaussian = 0;
           double tri_err_gaussian = 0;


            for (auto m : M)
                for (auto n : N)
                    for (int q = 0; q < 2; q++)
                    {
                        //! [curvature-measures-SurfaceMesh]
                        SM smesh;
                        double exp_H_min = 0.0;
                        double exp_H_max = 0.0;
                        double exp_G_min = 0.0;
                        double exp_G_max = 0.0;
                        if (input == "torus")
                        {
                            const double big_radius = 3.0;
                            const double small_radius = 1.0;
                            smesh = SMH::makeTorus(big_radius, small_radius,
                                RealPoint{ 0.0, 0.0, 0.0 }, m, n, 0,
                                SMH::NormalsType::VERTEX_NORMALS);
                            exp_H_min = (0.5 / (small_radius - big_radius) + 0.5 / small_radius);
                            exp_H_max = (0.5 / (big_radius + small_radius) + 0.5 / small_radius);
                            exp_G_min = (1.0 / (small_radius - big_radius) * 1.0 / small_radius);
                            exp_G_max = (1.0 / (big_radius + small_radius) * 1.0 / small_radius);
                        }
                        else if (input == "sphere")
                        {
                            const double radius = 2.0;
                            smesh = SMH::makeSphere(radius, RealPoint{ 0.0, 0.0, 0.0 }, m, n,
                                SMH::NormalsType::VERTEX_NORMALS);
                            exp_H_min = 1.0 / radius;
                            exp_H_max = 1.0 / radius;
                            exp_G_min = 1.0 / (radius * radius);
                            exp_G_max = 1.0 / (radius * radius);
                        }
                        else if (input == "lantern")
                        {
                            const double radius = 2.0;
                            smesh = SMH::makeLantern(radius, 1.0, RealPoint{ 0.0, 0.0, 0.0 }, m, n,
                                SMH::NormalsType::VERTEX_NORMALS);
                            exp_H_min = 0.5 / radius;
                            exp_H_max = 0.5 / radius;
                            exp_G_min = 0.0;
                            exp_G_max = 0.0;
                        }
                        //! [curvature-measures-SurfaceMesh]

                        //! [curvature-measures-CNC]
                        // builds a CorrectedNormalCurrentComputer object onto the torus/lantern/sphere mesh
                        CNC cnc(smesh);
                        // computes area, mean and Gaussian curvature measures
                        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
                        auto mu0 = cnc.computeMu0(q);
                        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
                        ((q) ? quad_time_0 : tri_time_0) += std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();

                        begin = std::chrono::steady_clock::now();
                        auto mu1 = cnc.computeMu1(q);
                        end = std::chrono::steady_clock::now();
                        ((q) ? quad_time_1 : tri_time_1) += std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();

                        begin = std::chrono::steady_clock::now();
                        auto mu2 = cnc.computeMu2(q);
                        end = std::chrono::steady_clock::now();
                        ((q) ? quad_time_2 : tri_time_2) += std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
                        //! [curvature-measures-CNC]

                        //! [curvature-measures-estimations]
                        // estimates mean (H) and Gaussian (G) curvatures by measure normalization.
                        std::vector< double > H(smesh.nbFaces());
                        std::vector< double > G(smesh.nbFaces());
                        for (auto f = 0; f < smesh.nbFaces(); ++f)
                        {
                            const auto b = smesh.faceCentroid(f);
                            const auto area = mu0.measure(b, r, f);
                            H[f] = cnc.meanCurvature(area, mu1.measure(b, r, f));
                            G[f] = cnc.GaussianCurvature(area, mu2.measure(b, r, f));
                        }
                        //! [curvature-measures-estimations]

                        //! [curvature-measures-check]
                        auto H_min_max = std::minmax_element(H.cbegin(), H.cend());
                        auto G_min_max = std::minmax_element(G.cbegin(), G.cend());
                        
                        ((q) ? quad_err_mean : tri_err_mean) += ( abs(exp_H_min - (*H_min_max.first)) + abs(exp_H_max - (*H_min_max.second)));
                        ((q) ? quad_err_gaussian : tri_err_gaussian) += ( abs(exp_G_min - (*G_min_max.first)) + abs(exp_G_max - (*G_min_max.second)));
  
                        //! [curvature-measures-check]

                        // [curvature-measures-output]
                        /*typedef SurfaceMeshWriter< RealPoint, RealVector > SMW;
                        const auto colormapH = makeQuantifiedColorMap(makeColorMap(-0.625, 0.625));
                        const auto colormapG = makeQuantifiedColorMap(makeColorMap(-0.625, 0.625));
                        auto colorsH = SMW::Colors(smesh.nbFaces());
                        auto colorsG = SMW::Colors(smesh.nbFaces());
                        for (auto i = 0; i < smesh.nbFaces(); i++)
                        {
                            colorsH[i] = colormapH(H[i]);
                            colorsG[i] = colormapG(G[i]);
                        }
                        SMW::writeOBJ("example-cnc-H"+input, smesh, colorsH);
                        SMW::writeOBJ("example-cnc-G"+input, smesh, colorsG);*/
                    }
                    std::cout << input << ": using m{10,20,50,75,100} and n{{10,20,50,75,100}" << "|| R:" << r << std::endl
                        << "evaluating measures on quads by\n"
                        << "(Barycentric Triangulation + Triangle Formulas) Vs. Bilinear Interpolat quad Formulas\n\n"

                        << "Sum of Errors:\n"
                        << "Mean: " << tri_err_mean << " Vs." << quad_err_mean << std::endl
                        << "Gaussian: " << tri_err_gaussian << " Vs." << quad_err_gaussian << std::endl << std::endl
                        << "___________________________________________________________\n\n";

        }
        std::cout << "Sum of run times for Mu_i:\n"
            << "mu0 tri: " << tri_time_0 << " Vs. quad:" << quad_time_0 << std::endl
            << "mu1 tri: " << tri_time_1 << " Vs. quad:" << quad_time_1 << std::endl
            << "mu2 tri: " << tri_time_2 << " Vs. quad:" << quad_time_2 << std::endl;
    system("pause");
}
