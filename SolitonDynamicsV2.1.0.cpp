//
//                                      SolitonDynamicsV2.0.cpp

//                              Created by Islam Albazlamit on 27/03/2025.
//           Original Code Created in Fortran by Prof. Usama Al Khawaja and Dr. Laila Al Sakkaf
//                       For more info please contact: islam.albazlamit@outlook.com

//i \psi_t(x,t) + a_1 \psi_{xx}(x,t) + a_2 |\psi(x,t)|^2 \psi(x,t) + a_3 |\psi(x,t)|^4 \psi(x,t) + V(x) \psi(x,t) = 0

// Solution used in this code is DOI: https://doi.org/10.1103/PhysRevE.107.014202

 #include <iostream>
 #include <vector>
 #include <cmath>
 #include <array>
 #include <complex>
 //#include <boost/multiprecision/cpp_dec_float.hpp>
 //#include <boost/multiprecision/complex_adaptor.hpp>
 //#include <boost/math/special_functions/pow.hpp>
 #include <chrono>
 #include <fstream>
 #include <cstdio>


 using namespace std;
 using namespace std::chrono;
 //using namespace boost::multiprecision;
 //using float_50 = cpp_dec_float_50;
 //using complex_50 = complex<float_50>;

int main(){
    
    int L = 400, n = 1000 /* # Segments in x */, np = 4, /* Skip in x grid */ npic = 200;
    int nn = 4, beta = 1 /* Sol Dir Right*/, tend = 100, numsoli,  nsteps;
    long double a0[n], b0[n], a1[n], b1[n],a2[n], b2[n], a3[n], b3[n], a4[n], b4[n], f[n], g[n], v[n];
    long double fp[n-8], gp[n-8];
    long double dx, g1, g2, g3, g30, x, x0, dt, norm, ke,ie;
    long double pe, Gamma, u0, V00, v0, /*v02(0),*/ t, u02, tot;
    long double A1, lambda1, lambda2;
    
    
    ofstream denout("density.dat");
    ofstream keout("ke.dat");
    ofstream peout("pe.dat");
    ofstream normout("norm.dat");
    ofstream ieout("ie.dat");
    ofstream animout("animation.dat");
    ofstream total("total.dat");
    ofstream vout("vout.dat");
    g1 = 0.5; g2=1; norm = 0;ke =0; ie=0;
    pe=0; Gamma = 0.000000000001;
    u0 = 0.5; V00 = 8; A1 =1;
    dx = L / (n - 1.0);
    dt = 0.5 * powl(dx, 2);             // discretization in t
    nsteps= (int)(tend/dt);
    g30 = 3.0 * powl(g2, 2) / (16.0 * u0);
    u0 = powl(g2, 2) * (nn + 1) / (g30 * powl((nn+2), 2));
    u02 = u0;
    g3 =  Gamma * g30;
    v0 =  -0.25 * 8;
    lambda1 = sqrtl(1+Gamma);
    lambda2 = sqrtl((2* u0)/(g2* lambda1));
    t = 0;
    
    auto start = steady_clock::now();
    for(int i = 0; i<= n ; i++){
        a0[i] = 0;
        b0[i] = 0;
        a1[i] = 0;
        b1[i] = 0;
        a2[i] = 0;
        b2[i] = 0;
        a3[i] = 0;
        b3[i] = 0;
        a4[i] = 0;
        b4[i] = 0;
        f[i] = 0;
        g[i] = 0;
        v[i] = 0;
    }
    
    for(int kkk=1; kkk <=1; kkk++){
        x0 = 30.0 * beta + 0.125 * kkk;
        cout << "x0 = " << x0 << endl;
        for(int iv = 1; iv <= 1; iv++){
            for(int jj = 1; jj <= 1; jj++){
                for(int i =1/*324*/; i <= n; i++){
                    x = (i - 1) * dx - L / 2.0;
                    numsoli = 0;        // 0 one soliton and 1 for two
                    f[i] = cosl((v0 / (2*g1)) * (x-x0) - ((v0*v0)/(4*g1))*t) * lambda2/(sqrtl(((1-lambda1)/2*lambda1) + powl(coshl((sqrtl(u0/g1) *(-(t * v0) + x - x0) )),2))) + numsoli *(cosl((v0 / (2*g1)) * (x-x0) - ((v0*v0)/(4*g1))*t) * lambda2/(sqrtl(((1-lambda1)/2*lambda1) + powl(coshl((sqrtl(u0/g1) *(-(t * v0) + x - x0) )),2))));
                    
                    g[i] = sinl((v0 / (2*g1)) * (x-x0) - ((v0*v0)/(4*g1))*t)  *lambda2/(sqrtl(((1-lambda1)/2*lambda1) + powl(coshl((sqrtl(u0/g1) *(-(t * v0) + x - x0) )),2)))+numsoli * (sinl((v0 / (2*g1)) * (x-x0) - ((v0*v0)/(4*g1))*t) * lambda2/(sqrtl(((1-lambda1)/2*lambda1) + powl(coshl((sqrtl(u0/g1) *(-(t * v0) + x - x0) )),2))));
                    // Test
                    //                    cout << setprecision(16);
                    //                    cout << i << " " << f[i] << " " << g[i] << endl;
                }
                
                for(int i =1; i <= n; i++){
                    x = (i - 1) * dx - L / 2.0;
                    v[i] = - 0* (V00) / powl((coshl( sqrtl(fabsl(V00))* (x))), 2);
                    
                    //                  To increase # of digits activate the following lines/add the # needed
                    //                    vout << setprecision(50);
                    vout << x << " " << v[i] << '\n';
                }
                
                
                for(int istep = 1; istep <= nsteps; istep ++){
                    t= istep * dt;
                        if ((fmod(istep, nsteps / npic) == 0) || (istep == 1)) {
                            for(int i = 1; i <= n ; i+=np){
                                x = (i-1) * dx - L / 2.0;
                                
                                //                  To increase # of digits activate the following lines/add the # needed
                                
                                denout << setprecision(16);
                                animout << setprecision(16);
                                denout << x << " " << t << " " << sqrtl(powl(f[i], 2.0) + powl(g[i], 2.0)) << '\n';
                                animout << x << " " << sqrtl(powl(f[i], 2.0) + powl(g[i], 2.0)) << '\n';
                                }
                            }
                    
                        
                        // 9-Points stencil //
                        
                        for(int i = 5; i <= n-4 ; i++){
                            fp[i] = (3.0 * f[i-4] - 32.0 * f[i-3] + 168.0 * f[i-2] - 672.0 * f[i-1]
                                     - 3.0 * f[i+4] + 32.0 * f[i+3] - 168.0 * f[i+2] + 672.0 * f[i+1])
                            / (840.0 * dx);
                            
                            gp[i]= (3.0 * g[i-4] - 32.0 * g[i-3] + 168.0 * g[i-2] - 672.0 * g[i-1]
                                    - 3.0 * g[i+4] + 32.0 * g[i+3] - 168.0 * g[i+2] + 672.0 * g[i+1])
                            / (840.0 * dx);
                        }
                        
                        // Kinetic Energy //
                        ke =0; pe=0; norm =0; ie =0;
                        for(int i =5; i <= n-4; i++){
                            ke = ke + dx * g1 * (fp[i] * fp[i] + gp[i] * gp[i]);
                        }
                        
                        // Potential Energy //
                        for(int i =1; i <= n-4; i++){
                            pe = pe + v[i] * (f[i] * f[i] + g[i] * g[i]);
                        }
                        
                        // Interaction Energy //
                        for(int i =1; i <= n-4; i++){
                            ie = ie + dx * ((-1.0 / 2.0) * g2 * (f[i] * f[i] + g[i] * g[i]) * (f[i] * f[i] + g[i] * g[i])
                                            + (-1.0 / 3.0) * g3 * (f[i] * f[i] + g[i] * g[i]) * (f[i] * f[i] + g[i] * g[i]) * (f[i] * f[i] + g[i] * g[i]));
                        }
                        
                        // Total Norm //
                        for(int i =5; i <= n-4; i++){
                            norm = norm + dx * (f[i] * f[i] + g[i] * g[i]);
                        }
                        
                        tot = ke + pe + ie;
                        //                  To increase # of digits activate the following lines/add the # needed
                        
                        //                        keout << setprecision(50);
                        //                        peout << setprecision(50);
                        //                        ieout << setprecision(50);
                        //                        normout << setprecision(50);
                        //                        total << setprecision(50);
                        keout << t << " " << ke<< '\n';
                        peout << t << " " << pe << '\n';
                        ieout << t << " " << ie<< '\n';
                        normout << t << " " << norm<< '\n';
                        total << t << " " << tot << '\n';
                        
                        for(int i = 2; i<= n-1 ; i++){
                            a0[i] = f[i];
                            b0[i] = g[i];
                        }
                        
                        for(int s = 12; s<= n-11; s++){
                            a1[s] =  -(g3 * a0[s] * a0[s] * a0[s] * a0[s]  * b0[s]) - g2 * b0[s] * b0[s] * b0[s] - g3 * b0[s] * b0[s] * b0[s] * b0[s] * b0[s]
                            - a0[s] * a0[s]  * b0[s] * (g2 + 2.0 * g3 * b0[s] * b0[s] )
                            - (g1 * (75600.0 * b0[-11 + s] - 10053996959110.0 * b0[s] -1331.0 *
                                     (1512.0 * b0[-10 + s] - 19600.0 * b0[-9 + s] + 165375.0 * b0[-8 + s]
                                      - 1026000.0 * b0[-7 + s] + 5027400.0 * b0[-6 + s] - 20511792.0 * b0[-5 + s]
                                      + 73256400.0 * b0[-4 + s] - 244188000.0 * b0[-3 + s] + 854658000.0 * b0[-2 + s]
                                      - 4444221600.0 * b0[-1 + s] - 4444221600.0 * b0[1 + s] + 854658000.0 * b0[2 + s]
                                      - 244188000.0 * b0[3 + s] + 73256400.0 * b0[4 + s] - 20511792.0 * b0[5 + s]
                                      + 25.0 * (201096.0 * b0[6 + s] - 41040.0 * b0[7 + s] + 6615.0 * b0[8 + s]
                                                - 784.0 * b0[9 + s]) + 1512.0 * b0[10 + s]) + 75600.0 * b0[11 + s]))
                            / ((3226504881600.0 * dx * dx) - b0[s] * v[s]);
                            
                            
                        }
                        for(int s = 12; s<= n-11; s++){
                            
                            b1[s] = (g3 * a0[s] * a0[s] * a0[s] * a0[s] * a0[s] ) + (g1 * (75600.0 * a0[-11 + s] - 10053996959110.0 * a0[s]
                                                                                           - 1331.0 * (1512.0 * a0[-10 + s] - 19600.0 * a0[-9 + s] + 165375.0 * a0[-8 + s]
                                                                                                       - 1026000.0 * a0[-7 + s] + 5027400.0 * a0[-6 + s] - 20511792.0 * a0[-5 + s]
                                                                                                       + 73256400.0 * a0[-4 + s] - 244188000.0 * a0[-3 + s] + 854658000.0 * a0[-2 + s]
                                                                                                       - 4444221600.0 * a0[-1 + s] - 4444221600.0 * a0[1 + s] + 854658000.0 * a0[2 + s]
                                                                                                       - 244188000.0 * a0[3 + s] + 73256400.0 * a0[4 + s] - 20511792.0 * a0[5 + s]
                                                                                                       + 25.0 * (201096.0 * a0[6 + s] - 41040.0 * a0[7 + s] + 6615.0 * a0[8 + s]
                                                                                                                 - 784.0 * a0[9 + s]) + 1512.0 * a0[10 + s]) + 75600.0 * a0[11 + s]))
                            / (3226504881600.0 * dx * dx) + a0[s] * a0[s] * a0[s] * (g2 + 2.0 * g3 * b0[s] * b0[s] )
                            + a0[s] * (g2 * b0[s] * b0[s]  + g3 * b0[s] * b0[s] * b0[s] * b0[s]  + v[s]);
                            
                        }
                        for(int s = 12; s<= n-11; s++){
                            a2[s] = (-4.0 * g3 * a0[s] * a0[s] * a0[s]  * a1[s] * b0[s] - 2.0 * a0[s] * a1[s] * b0[s] * (g2 + 2.0 * g3 * b0[s] * b0[s] ) - g3 * a0[s] * a0[s] * a0[s] * a0[s]  * b1[s] - 3.0 * g2 * b0[s] * b0[s]  * b1[s]
                                     - 5.0 * g3 * b0[s] * b0[s] * b0[s] * b0[s]  * b1[s] - a0[s] * a0[s]  * (g2 + 6.0 * g3 * b0[s] * b0[s] )
                                     * b1[s] - (g1 * (75600.0 * b1[-11 + s] - 10053996959110.0 * b1[s] -1331.0 * (1512.0
                                                                                                                  * b1[-10 + s] - 19600.0 * b1[-9 + s] + 165375.0 * b1[-8 + s] - 1026000.0 * b1[-7 + s]
                                                                                                                  + 5027400.0 * b1[-6 + s] - 20511792.0 * b1[-5 + s] + 73256400.0 * b1[-4 + s] - 244188000.0
                                                                                                                  * b1[-3 + s] + 854658000.0 * b1[-2 + s] - 4444221600.0 * b1[-1 + s] - 4444221600.0
                                                                                                                  * b1[1 + s] + 854658000.0 * b1[2 + s] - 244188000.0 * b1[3 + s] + 73256400.0 * b1[4 + s]
                                                                                                                  - 20511792.0 * b1[5 + s] + 25.0 * (201096.0 * b1[6 + s] - 41040.0 * b1[7 + s] + 6615.0
                                                                                                                                                     * b1[8 + s] - 784.0 * b1[9 + s]) + 1512.0 * b1[10 + s]) + 75600.0 * b1[11 + s]))
                                     / (3226504881600.0 * dx * dx) - b1[s] * v[s]) / 2.0;
                            
                        }
                        for(int s = 12; s<= n-11; s++){
                            b2[s] = (5.0 * g3 * a0[s] * a0[s] * a0[s] * a0[s]  * a1[s] + (g1 * (75600.0 * a1[-11 + s] - 10053996959110.0
                                                                                                * a1[s] - 1331.0 * (1512.0 * a1[-10 + s] - 19600.0 * a1[-9 + s] + 165375.0 * a1[-8 + s]
                                                                                                                    - 1026000.0 * a1[-7 + s] + 5027400.0 * a1[-6 + s] - 20511792.0 * a1[-5 + s] + 73256400.0
                                                                                                                    * a1[-4 + s] - 244188000.0 * a1[-3 + s] + 854658000.0 * a1[-2 + s] - 4444221600.0
                                                                                                                    * a1[-1 + s] - 4444221600.0 * a1[1 + s] + 854658000.0 * a1[2 + s] - 244188000.0 * a1[3 + s]
                                                                                                                    + 73256400.0 * a1[4 + s] - 20511792.0 * a1[5 + s] + 25.0 * (201096.0 * a1[6 + s]
                                                                                                                                                                                - 41040.0 * a1[7 + s] + 6615.0 * a1[8 + s] - 784.0 * a1[9 + s]) + 1512.0 * a1[10 + s])
                                                                                                + 75600.0 * a1[11 + s])) / (3226504881600.0 * dx * dx) + 3.0 * a0[s] * a0[s]  * a1[s]
                                     * (g2 + 2.0 * g3 * b0[s] * b0[s] ) + 4.0 * g3 * a0[s] * a0[s] * a0[s]  * b0[s] * b1[s] + 2.0 * a0[s]
                                     * b0[s] * (g2 + 2.0 * g3 * b0[s] * b0[s] ) * b1[s] + a1[s] * (g2 * b0[s] * b0[s]  + g3
                                                                                                   * b0[s] * b0[s] * b0[s] * b0[s]  + v[s])) / 2.0;
                        }
                        for(int s = 12; s<= n-11; s++){
                            a3[s] = (-4.0 * g3 * a0[s] * a0[s] * a0[s]  * a2[s] * b0[s] - 2.0 * a0[s] * a2[s] * b0[s] * (g2 + 2.0 * g3
                                                                                                                         * b0[s] * b0[s] ) - a1[s] * a1[s]  * b0[s] * (g2 + 6.0 * g3 * a0[s] * a0[s]  + 2.0 * g3
                                                                                                                                                                       * b0[s] * b0[s] ) - 2.0 * a0[s] * a1[s] * (g2 + 2.0 * g3 * a0[s] * a0[s]  + 6.0 * g3
                                                                                                                                                                                                                  * b0[s] * b0[s] ) * b1[s] - 3.0 * g2 * b0[s] * b1[s] * b1[s]  - 10.0 * g3 * b0[s] * b0[s] * b0[s]
                                     * b1[s] * b1[s]  - g3 * a0[s] * a0[s] * a0[s] * a0[s]  * b2[s] - 3.0 * g2 * b0[s] * b0[s]  * b2[s]
                                     - 5.0 * g3 * b0[s] * b0[s] * b0[s] * b0[s]  * b2[s] - a0[s] * a0[s]  * (6.0 * g3 * b0[s] * b1[s] * b1[s]
                                                                                                             + g2 * b2[s] + 6.0 * g3 * b0[s] * b0[s]  * b2[s]) - (g1 * (75600.0 * b2[-11 + s]
                                                                                                                                                                        - 10053996959110.0 * b2[s] - 1331.0 * (1512.0 * b2[-10 + s] - 19600.0 * b2[-9 + s]
                                                                                                                                                                                                               + 165375.0 * b2[-8 + s] - 1026000.0 * b2[-7 + s] + 5027400.0 * b2[-6 + s] - 20511792.0
                                                                                                                                                                                                               * b2[-5 + s] + 73256400.0 * b2[-4 + s] - 244188000.0 * b2[-3 + s] + 854658000.0 * b2[-2 + s]
                                                                                                                                                                                                               - 4444221600.0 * b2[-1 + s] - 4444221600.0 * b2[1 + s] + 854658000.0 * b2[2 + s]
                                                                                                                                                                                                               - 244188000.0 * b2[3 + s] + 73256400.0 * b2[4 + s] - 20511792.0 * b2[5 + s] + 25.0
                                                                                                                                                                                                               * (201096.0 * b2[6 + s] - 41040.0 * b2[7 + s] + 6615.0 * b2[8 + s] - 784.0 * b2[9 + s])
                                                                                                                                                                                                               + 1512.0 * b2[10 + s]) + 75600.0 * b2[11 + s]))/ (3226504881600.0 * dx * dx) - b2[s]
                                     * v[s]) / 3.0;
                            
                        }
                        for(int s = 12; s<= n-11; s++){
                            
                            b3[s] = (5.0 * g3 *  a0[s] * a0[s] * a0[s] * a0[s]  *  a2[s] + (g1 * (75600.0 *  a2[-11 + s] - 10053996959110.0
                                                                                                  * a2[s] - 1331.0 * (1512.0 *  a2[-10 + s] - 19600.0 *  a2[-9 + s] + 165375.0 *  a2[-8 + s]
                                                                                                                      - 1026000.0 * a2[-7 + s] + 5027400.0 * a2[-6 + s] - 20511792.0 *  a2[-5 + s] + 73256400.0
                                                                                                                      * a2[-4 + s] - 244188000.0 *  a2[-3 + s] + 854658000.0 * a2[-2 + s] - 4444221600.0
                                                                                                                      * a2[-1 + s] - 4444221600.0 *  a2[1 + s] + 854658000.0 * a2[2 + s] - 244188000.0
                                                                                                                      * a2[3 + s] + 73256400.0 * a2[4 + s] - 20511792.0 * a2[5 + s] + 25.0 * (201096.0
                                                                                                                                                                                              * a2[6 + s] - 41040.0 * a2[7 + s] + 6615.0 * a2[8 + s] - 784 *  a2[9 + s]) + 1512.0
                                                                                                                      * a2[10 + s]) + 75600.0 * a2[11 + s]))/ (3226504881600.0 * dx * dx) + 2.0 * g2 * a1[s]
                                     * b0[s]* b1[s] + 4.0 * g3 *  a1[s]* b0[s] * b0[s] * b0[s]  * b1[s] + 3.0 * a0[s] * a0[s]  * (a2[s]
                                                                                                                                  * (g2 + 2* g3 * b0[s] * b0[s] ) + 4.0 * g3 * a1[s] * b0[s]* b1[s]) + 2.0 * g3 * a0[s] * a0[s] * a0[s]
                                     * (5.0 *  powl(a1[s], 2) +  b1[s] * b1[s]  + 2.0 * b0[s]* b2[s]) + a0[s] * (3.0 * powl(a1[s], 2)
                                                                                                                 * (g2 + 2.0 * g3 * b0[s] * b0[s] ) + (g2 + 6 * g3 * b0[s] * b0[s] ) * powl(b1[s],2) + 2.0
                                                                                                                 * b0[s]*(g2 + 2* g3 * b0[s] * b0[s] )* b2[s]) + a2[s] * (g2 * b0[s] * b0[s]  + g3
                                                                                                                                                                          * b0[s] * b0[s] * b0[s] * b0[s]  + v[s])) / 3.0;
                            
                        }
                        for(int s = 12; s<= n-11; s++){
                            a4[s] = (-2.0 * a1[s] * a2[s] * b0[s] * (g2 + 2.0 * g3 * b0[s] * b0[s] ) - powl(a1[s], 2) * (g2 + 6.0
                                                                                                                         * g3 * b0[s] * b0[s] ) * b1[s] - g2 * b1[s] * b1[s] * b1[s]  - 10.0 * g3 * b0[s] * b0[s]  * b1[s] * b1[s] * b1[s]
                                     - 6.0 * g2 * b0[s] * b1[s]* b2[s] - 20.0 * g3 * b0[s] * b0[s] * b0[s]  * b1[s]* b2[s] - 4.0 * g3
                                     * a0[s] * a0[s] * a0[s] * (a3[s] * b0[s] + a2[s] * b1[s] + a1[s] * b2[s]) - 2.0 * a0[s] * (2.0 * g3
                                                                                                                                * powl(a1[s], 3) * b0[s] + a3[s] * b0[s] * (g2 + 2.0 * g3 * b0[s] * b0[s] ) + a2[s] * (g2 + 6.0
                                                                                                                                                                                                                       * g3 * b0[s] * b0[s] )* b1[s] + a1[s] * (6.0 * g3 * b0[s] * b1[s] * b1[s]  + g2 * b2[s] + 6.0
                                                                                                                                                                                                                                                                * g3 * b0[s] * b0[s]  * b2[s])) - g3 * a0[s] * a0[s] * a0[s] * a0[s]  * b3[s] - 3 * g2 * b0[s] * b0[s]  * b3[s]
                                     - 5.0 * g3 * b0[s] * b0[s] * b0[s] * b0[s]  * b3[s] - a0[s] * a0[s]  * (12 * g3 *  a1[s] * a2[s] * b0[s]
                                                                                                             + 6.0 * g3 * powl(a1[s], 2) * b1[s] + 2.0 * g3 * b1[s] * b1[s] * b1[s]  + 12.0 * g3 * b0[s] * b1[s]
                                                                                                             * b2[s] + (g2 + 6.0 * g3 * b0[s] * b0[s] ) * b3[s]) - (g1 * (75600.0 * b3[-11 + s]
                                                                                                                                                                          - 10053996959110.0 * b3[s] - 1331.0 * (1512.0 * b3[-10 + s] - 19600.0 *b3[-9 + s] + 165375.0
                                                                                                                                                                                                                 * b3[-8 + s] - 1026000.0 * b3[-7 + s] + 5027400.0 * b3[-6 + s] - 20511792.0
                                                                                                                                                                                                                 * b3[-5 + s] + 73256400.0 *b3[-4 + s] - 244188000.0 * b3[-3 + s] + 854658000.0 * b3[-2 + s]
                                                                                                                                                                                                                 - 4444221600.0 * b3[-1 + s] - 4444221600.0 * b3[1 + s] + 854658000.0 * b3[2 + s] - 244188000.0
                                                                                                                                                                                                                 * b3[3 + s] + 73256400.0 *b3[4 + s] - 20511792.0 *b3[5 + s] + 25.0 *(201096.0
                                                                                                                                                                                                                                                                                      * b3[6 + s] - 41040.0 *b3[7 + s] + 6615.0 *b3[8 + s] - 784.0 * b3[9 + s]) + 1512.0
                                                                                                                                                                                                                 * b3[10 + s]) + 75600.0 * b3[11 + s]))/ (3226504881600.0 * powl(dx, 2)) -  b3[s]
                                     * v[s])/4.0;
                            
                        }
                        for(int s = 12; s<= n-11; s++){
                            
                            b4[s] = (5.0 * g3 * a0[s] * a0[s] * a0[s] * a0[s]  * a3[s] + (g1 * (75600.0 * a3[-11 + s] - 10053996959110.0 * a3[s]
                                                                                                - 1331.0 * (1512.0 * a3[-10 + s] - 19600.0 * a3[-9 + s] + 165375.0 * a3[-8 + s] - 1026000.0
                                                                                                            * a3[-7 + s] + 5027400.0 * a3[-6 + s] - 20511792.0 * a3[-5 + s] + 73256400.0 * a3[-4 + s]
                                                                                                            - 244188000.0 * a3[-3 + s] + 854658000.0 * a3[-2 + s] - 4444221600.0 * a3[-1 + s]
                                                                                                            - 4444221600.0 * a3[1 + s] + 854658000.0 * a3[2 + s] - 244188000.0 * a3[3 + s] + 73256400.0
                                                                                                            * a3[4 + s] - 20511792.0 * a3[5 + s] + 25.0 * (201096.0 * a3[6 + s] - 41040.0 * a3[7 + s]
                                                                                                                                                           + 6615.0 * a3[8 + s] - 784.0 * a3[9 + s]) + 1512.0 * a3[10 + s]) + 75600.0 * a3[11 + s]))
                                     / (3226504881600.0 * powl(dx, 2)) + g2 * a3[s] * b0[s] * b0[s]  + g3 * a3[s]
                                     * b0[s] * b0[s] * b0[s] * b0[s]  + powl(a1[s], 3) * (g2 + 10.0 * g3 * a0[s] * a0[s]  + 2 * g3 * b0[s] * b0[s] )
                                     + 12.0 * g3 *  a0[s]*  powl(a1[s], 2) * b0[s] * b1[s] + 2.0 * g2 * a2[s]* b0[s]* b1[s]
                                     + 4.0 * g3 *  a2[s]* b0[s] * b0[s] * b0[s]  * b1[s] + 3.0 *  a0[s] * a0[s] *( a3[s]*(g2 + 2* g3
                                                                                                                          * b0[s] * b0[s] ) + 4.0 * g3 *  a2[s]* b0[s]* b1[s]) + a1[s]*(20.0 * g3 * a0[s] * a0[s] * a0[s]
                                                                                                                                                                                        * a2[s] + 6.0 *  a0[s]*  a2[s] * (g2 + 2 * g3 * b0[s] * b0[s] ) + (g2 + 6 * g3 * b0[s] * b0[s] )
                                                                                                                                                                                        * b1[s] * b1[s]  + 2.0 * b0[s] * (g2 + 2 * g3 * b0[s] * b0[s] ) * b2[s] + 6.0 * g3
                                                                                                                                                                                        * a0[s] * a0[s]  *( b1[s] * b1[s]  + 2 * b0[s] * b2[s])) + 4 * g3 *  a0[s] * a0[s] * a0[s]
                                     * (b1[s]* b2[s] + b0[s]* b3[s]) + 2.0 *  a0[s] * (g2* b1[s]* b2[s] + 6.0 * g3 * b0[s] * b0[s]
                                                                                       * b1[s]* b2[s] + 2.0 * g3 * b0[s] * b0[s] * b0[s]  * b3[s] + b0[s]*(2.0 * g3 * b1[s] * b1[s] * b1[s]
                                                                                                                                                           + g2 * b3[s])) + a3[s]* v[s]) / 4.0;
                            
                        }
                        
                        // Updating System //
                        for(int i = 1; i <= n; i++){
                            f[i] = a0[i] + a1[i] * dt + a2[i] * powl(dt, 2) + a3[i] * powl(dt, 3)
                            + a4[i] * powl(dt, 4);
                            g[i] = b0[i] + b1[i] * dt + b2[i] * powl(dt, 2) + b3[i] * powl(dt, 3)
                            + b4[i] * powl(dt, 4);
                            
                            // Exit when soliton reaches boundaries //
                        }
                        
                        //                    t= istep * dt;
                    }
                    cout << "End of istep" << '\n';
                }
            cout << "jj" << '\n';
        }
                
        auto end = chrono::steady_clock::now();
        auto diff = end - start;
        cout << (duration <double, milli>(diff).count())/ 1000 << "s" << '\n';
        cout << (duration <double, milli>(diff).count())/ 60000 << "mins" << '\n';
        
        denout.close();
        keout.close();
        peout.close();
        ieout.close();
        normout.close();
        denout.close();
        total.close();

        return 0;
    }
}

