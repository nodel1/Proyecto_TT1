// $Source$
//--------------------------------------------------
// DEInteg
//--------------------------------------------------
// Proyecto_TT1: Numerical Integration for ODEs
//
// Copyright (c) 2025, Meysam Mahooti
//
// Created: 2025/05/23
//
/** @file DEInteg.cpp
 *  @brief Numerical integration for ordinary differential equations using Shampine & Gordon's variable order, variable stepsize multistep method
 *
 *  @author NOEL DEL RIO GONZALEZ
 *  @bug No known bugs
 */

#include <cmath>
#include <limits>
#include <stdexcept>
#include <iostream>
#include "..\include\DEInteg.hpp"



Matrix DEInteg(Matrix (*f)(double, Matrix&), double t, double tout, double relerr, double abserr, int n_eqn, Matrix& y) {
    // Asegurar que y sea una matriz columna
    if (y.n_row < y.n_column) {
        y = transpose(y);
    }
	
	cout << "primero del dei" << endl;

    // Constantes de precisión de máquina
    const double twou = 2.0 * std::numeric_limits<double>::epsilon();     //no viene en satconst pero lo mismo lo suyo era meterel epsilon ahi revisar luego
    const double fouru = 4.0 * std::numeric_limits<double>::epsilon();

    // Estados del integrador
    struct {
        int DE_INIT = 1;      
        int DE_DONE = 2;      
        int DE_BADACC = 3;    
        int DE_NUMSTEPS = 4;  
        int DE_STIFF = 5;     
        int DE_INVPARAM = 6;  
    } DE_STATE;

    int State_ = DE_STATE.DE_INIT;
    bool PermitTOUT = true; 
    double told = 0.0;



	cout << "segundo del dei" << endl;


    // Inicializar vectores constantes
    Matrix two(14, 1);
    two(1,1) = 1.0; two(2,1) = 2.0; two(3,1) = 4.0; two(4,1) = 8.0; two(5,1) = 16.0;
    two(6,1) = 32.0; two(7,1) = 64.0; two(8,1) = 128.0; two(9,1) = 256.0;
    two(10,1) = 512.0; two(11,1) = 1024.0; two(12,1) = 2048.0; two(13,1) = 4096.0; two(14,1) = 8192.0;

    Matrix gstr(14, 1);
    gstr(1,1) = 1.0; gstr(2,1) = 0.5; gstr(3,1) = 0.0833; gstr(4,1) = 0.0417; gstr(5,1) = 0.0264;
    gstr(6,1) = 0.0188; gstr(7,1) = 0.0143; gstr(8,1) = 0.0114; gstr(9,1) = 0.00936; gstr(10,1) = 0.00789;
    gstr(11,1) = 0.00679; gstr(12,1) = 0.00592; gstr(13,1) = 0.00524; gstr(14,1) = 0.00468;

    // Inicializar matrices de trabajo
    Matrix yy(n_eqn, 1);
    Matrix wt(n_eqn, 1);
    Matrix p(n_eqn, 1);
    Matrix yp(n_eqn, 1);
    Matrix phi(n_eqn, 17);
    Matrix g(14, 1);
    Matrix sig(14, 1);
    Matrix rho(14, 1);
    Matrix w(13, 1);
    Matrix alpha(13, 1);
    Matrix beta(13, 1);
    Matrix v(13, 1);
    Matrix psi_(13, 1);

    // Si t == tout, no se necesita integración
    if (t == tout) {
        return y;
    }

    // Verificar parámetros de entrada
    double epsilon = std::max(relerr, abserr);
    if (relerr < 0.0 || abserr < 0.0 || epsilon <= 0.0 || 
        State_ > DE_STATE.DE_INVPARAM || 
        (State_ != DE_STATE.DE_INIT && t != told)) {
        State_ = DE_STATE.DE_INVPARAM;
        throw std::invalid_argument("Parámetros inválidos en DEInteg");
    }
	
	
	
		cout << "tercero del dei" << endl;
		
		
		

    // Configurar intervalo de integración
    double del = tout - t;
    double absdel = std::abs(del);
    double tend = t + 100.0 * del;
    if (!PermitTOUT) {
        tend = tout;
    }

    int nostep = 0;
    int kle4 = 0;
    bool stiff = false;
    double releps = relerr / epsilon;
    double abseps = abserr / epsilon;
    bool start = true;
    double x = t;
    double delsgn = (del >= 0.0) ? 1.0 : -1.0;
    double h = sign_(std::max(fouru * std::abs(x), std::abs(tout - x)), tout - x);
    bool OldPermit = false;
    double hold = 0.0;
    int k = 1;
    int kold = 0;
    bool phase1 = true;
    bool nornd = true;
    int ns = 0;
    int ifail = 0;

    yy = y;
	
		cout << "cuarto del dei" << endl;
		
		

    while (true) { // Bucle principal de pasos
        // Interpolación si se pasa el punto de salida
        if (std::abs(x - t) >= absdel) {
            Matrix yout(n_eqn, 1);
            Matrix ypout(n_eqn, 1);
            g(2,1) = 1.0;                 //REVISAR ESTO por si es el uno o el dos de indice
            rho(2,1) = 1.0;
            double hi = tout - x;
            int ki = kold + 1;

            // Inicializar w para calcular g
            for (int i = 1; i <= ki; ++i) {
                w(i+1,1) = 1.0 / i;
            }
			
				cout << "quinto del dei" << endl;

            // Calcular g y rho
            double term = 0.0;
            for (int j = 2; j <= ki; j++) {            //REVISAR SI PONER j++ o ++j
                double psijm1 = psi_(j,1);
                double gamma = (hi + term) / psijm1;
                double eta = hi / psijm1;
                for (int i = 1; i <= ki + 1 - j; ++i) {
                    w(i+1,1) = gamma * w(i+1,1) - eta * w(i+2,1);
                }
                g(j+1,1) = w(2,1);
                rho(j+1,1) = gamma * rho(j,1);
                term = psijm1;
            }

            // Interpolar solución y derivada
            for (int j = 1; j <= ki; j++) {
                int i = ki + 1 - j;
                for (int l = 1; l <= n_eqn; ++l) {
                    yout(l,1) += g(i+1,1) * phi(l,i+1);
                    ypout(l,1) += rho(i+1,1) * phi(l,i+1);
                }
            }
			
				cout << "sexto del dei" << endl;
				
				
            yout = y + yout * hi;
            y = yout;
            State_ = DE_STATE.DE_DONE;
            t = tout;
            told = t;
            OldPermit = PermitTOUT;
            return y;
        }



	cout << "septimo del dei" << endl;
	
	
	
        // Extrapolación si está cerca de tout
		if (!PermitTOUT && std::abs(tout - x) < fouru * std::abs(x)) {
			h = tout - x;
			Matrix temp_yp = f(x, yy); // Matriz auxiliar
			yp = temp_yp;              // Asignar a yp
			y = yy + yp * h;
			State_ = DE_STATE.DE_DONE;
			t = tout;
			told = t;
			OldPermit = PermitTOUT;
			return y;
		}
		
			cout << "octavo del dei" << endl;

        // Limitar tamaño del paso y calcular pesos
        h = sign_(std::min(std::abs(h), std::abs(tend - x)), h);
        for (int l = 1; l <= n_eqn; l++) {
            wt(l,1) = releps * std::abs(yy(l,1)) + abseps;
            if (wt(l,1) == 0.0) {
                State_ = DE_STATE.DE_INVPARAM;
                throw std::runtime_error("Peso cero detectado en la estimación de error");
            }
        }
		
			cout << "noveno del dei" << endl;

        // Bloque 0: Verificar tamaño del paso y tolerancia
        bool crash = false;
        if (std::abs(h) < fouru * std::abs(x)) {
            h = sign_(fouru * std::abs(x), h);
            crash = true;
        }


	cout << "decinmo del dei" << endl;
	
	
        double p5eps = 0.5 * epsilon;
        g(2,1) = 1.0;
        g(3,1) = 0.5;
        sig(2,1) = 1.0;

        double round = 0.0;
        for (int l = 1; l <= n_eqn; l++) {
            round += (y(l,1) * y(l,1)) / (wt(l,1) * wt(l,1));
        }
        round = twou * std::sqrt(round);
        if (p5eps < round) {
            epsilon = 2.0 * round * (1.0 + fouru);
            crash = true;
        }

        if (crash) {
            State_ = DE_STATE.DE_BADACC;
            relerr = epsilon * releps;
            abserr = epsilon * abseps;
            y = yy;
            t = x;
            told = t;
            OldPermit = true;
            return y;
        }
		
		
			cout << "onceavo del dei" << endl;

        // Inicialización para el primer paso
        double absh;
        double hnew = 0.0;
        int knew = 0;
        if (start) {
			Matrix temp_yp = f(x, yy); // Matriz auxiliar
			yp = temp_yp;
            double sum = 0.0;
            for (int l = 1; l <= n_eqn; l++) {       //igual lo del ++ y revisar abajo
                phi(l,2) = yp(l,1);
                phi(l,3) = 0.0;
                sum += (yp(l,1) * yp(l,1)) / (wt(l,1) * wt(l,1));
            }
            sum = std::sqrt(sum);
            absh = std::abs(h);
            if (epsilon < 16.0 * sum * h * h) {
                absh = 0.25 * std::sqrt(epsilon / sum);
            }
            h = sign_(std::max(absh, fouru * std::abs(x)), h);
            hold = 0.0;
            kold = 0;
            start = false;
            nornd = (p5eps > 100.0 * round);
            if (!nornd) {
                for (int l = 1; l <= n_eqn; ++l) {
                    phi(l,16) = 0.0;
                }
            }
        }
		
			cout << "doceavo del dei" << endl;

        // Bucle para pasos exitosos
        int kp1, kp2, km1, km2;
        double erkm1 = 0.0, erkm2 = 0.0, erk = 0.0, erkp1 = 0.0;
        while (true) {
            // Bloque 1: Calcular coeficientes
            kp1 = k + 1;
            kp2 = k + 2;
            km1 = k - 1;
            km2 = k - 2;

            if (h != hold) {
                ns = 0;
            }
            if (ns <= kold) {
                ns++;
            }
            int nsp1 = ns + 1;

            if (k >= ns) {
                beta(ns+1,1) = 1.0;
                alpha(ns+1,1) = 1.0 / ns;
                double temp1 = h * ns;
                sig(nsp1+1,1) = 1.0;
                if (k >= nsp1) {
                    for (int i = nsp1; i <= k; ++i) {
                        int im1 = i - 1;
                        double temp2 = psi_(im1+1,1);
                        psi_(im1+1,1) = temp1;
                        beta(i+1,1) = beta(im1+1,1) * psi_(im1+1,1) / temp2;
                        temp1 = temp2 + h;
                        alpha(i+1,1) = h / temp1;
                        sig(i+2,1) = i * alpha(i+1,1) * sig(i+1,1);
                    }
                }
                psi_(k+1,1) = temp1;

                if (ns > 1) {
                    if (k > kold) {
                        double temp4 = k * kp1;
                        v(k+1,1) = 1.0 / temp4;
                        int nsm2 = ns - 2;
                        for (int j = 1; j <= nsm2; j++) {
                            int i = k - j;
                            v(i+1,1) = v(i+1,1) - alpha(j+2,1) * v(i+2,1);
                        }
                    }
                    int limit1 = kp1 - ns;
                    double temp5 = alpha(ns+1,1);
                    for (int iq = 1; iq <= limit1; iq++) {
                        v(iq+1,1) = v(iq+1,1) - temp5 * v(iq+2,1);
                        w(iq+1,1) = v(iq+1,1);
                    }
                    g(nsp1+1,1) = w(2,1);
                } else {
                    for (int iq = 1; iq <= k; iq++) {
                        double temp3 = iq * (iq + 1);
                        v(iq+1,1) = 1.0 / temp3;
                        w(iq+1,1) = v(iq+1,1);
                    }
                }

                int nsp2 = ns + 2;
                if (kp1 >= nsp2) {
                    for (int i = nsp2; i <= kp1; i++) {
                        int limit2 = kp2 - i;
                        double temp6 = alpha(i,1);
                        for (int iq = 1; iq <= limit2; iq++) {
                            w(iq+1,1) -= temp6 * w(iq+2,1);
                        }
                        g(i+1,1) = w(2,1);
                    }
                }
            }
			
				cout << "trece del dei" << endl;

            // Bloque 2: Predecir solución y estimar error
            for (int i = nsp1; i <= k; i++) {
                double temp1 = beta(i+1,1);
                for (int l = 1; l <= n_eqn; l++) {
                    phi(l,i+1) = temp1 * phi(l,i+1);
                }
            }

            for (int l = 1; l <= n_eqn; ++l) {
                phi(l,kp2+1) = phi(l,kp1+1);
                phi(l,kp1+1) = 0.0;
                p(l,1) = 0.0;
            }

            for (int j = 1; j <= k; ++j) {
                int i = kp1 - j;
                int ip1 = i + 1;
                double temp2 = g(i+1,1);
                for (int l = 1; l <= n_eqn; ++l) {
                    p(l,1) += temp2 * phi(l,i+1);
                    phi(l,i+1) += phi(l,ip1+1);
                }
            }

            if (nornd) {
                p = y + p * h;
            } else {
                for (int l = 1; l <= n_eqn; ++l) {
                    double tau = h * p(l,1) - phi(l,16);
                    p(l,1) = y(l,1) + tau;
                    phi(l,17) = (p(l,1) - y(l,1)) - tau;
                }
            }
			
				cout << "14 del dei" << endl;

            double xold = x;
            x += h;
            absh = std::abs(h);
            Matrix temp_yp = f(x, p); // Matriz auxiliar
			yp = temp_yp;

            // Estimar errores en órdenes k, k-1, k-2
            erkm2 = 0.0;
            erkm1 = 0.0;
            erk = 0.0;
            for (int l = 1; l <= n_eqn; ++l) {
                double temp3 = 1.0 / wt(l,1);
                double temp4 = yp(l,1) - phi(l,2);
                if (km2 > 0) {
                    erkm2 += std::pow((phi(l,km1+1) + temp4) * temp3, 2);
                }
                if (km2 >= 0) {
                    erkm1 += std::pow((phi(l,k+1) + temp4) * temp3, 2);
                }
                erk += std::pow(temp4 * temp3, 2);
            }

            if (km2 > 0) {
                erkm2 = absh * sig(km1+1,1) * gstr(km2+1,1) * std::sqrt(erkm2);
            }
            if (km2 >= 0) {
                erkm1 = absh * sig(k+1,1) * gstr(km1+1,1) * std::sqrt(erkm1);
            }

            double temp5 = absh * std::sqrt(erk);
            double err = temp5 * (g(k+1,1) - g(kp1+1,1));
            erk = temp5 * sig(kp1+1,1) * gstr(k+1,1);
            knew = k;

            if (km2 > 0 && std::max(erkm1, erkm2) <= erk) {
                knew = km1;
            }
            if (km2 == 0 && erkm1 <= 0.5 * erk) {
                knew = km1;
            }


			cout << "15 del dei" << endl;
	
	
            bool success = (err <= epsilon);

            // Bloque 3: Manejar paso fallido
            if (!success) {
                phase1 = false;
                x = xold;
                for (int i = 1; i <= k; i++) {
                    double temp1 = 1.0 / beta(i+1,1);
                    int ip1 = i + 1;
                    for (int l = 1; l <= n_eqn; ++l) {
                        phi(l,i+1) = temp1 * (phi(l,i+1) - phi(l,ip1+1));
                    }
                }

                if (k >= 2) {
                    for (int i = 2; i <= k; i++) {
                        psi_(i,1) = psi_(i+1,1) - h;
                    }
                }

                ifail++;
                double temp2 = 0.5;
                if (ifail > 3 && p5eps < 0.25 * erk) {
                    temp2 = std::sqrt(p5eps / erk);
                }
                if (ifail >= 3) {
                    knew = 1;
                }
                h *= temp2;
                k = knew;

                if (std::abs(h) < fouru * std::abs(x)) {
                    State_ = DE_STATE.DE_BADACC;
                    relerr = epsilon * releps;
                    abserr = epsilon * abseps;
                    y = yy;
                    t = x;
                    told = t;
                    OldPermit = true;
                    return y;
                }
            } else {
                break; // Paso exitoso, salir del bucle interno
            }
        }
		
		cout << "16 del dei" << endl;
			
			

        // Bloque 4: Corregir solución y actualizar diferencias
        kold = k;
        hold = h;

        double temp1 = h * g(kp1+1,1);
        for (int l = 1; l <= n_eqn; l++) {
            if (nornd) {
                y(l,1) = p(l,1) + temp1 * (yp(l,1) - phi(l,2));
            } else {
                double rho_val = temp1 * (yp(l,1) - phi(l,2)) - phi(l,17);
                y(l,1) = p(l,1) + rho_val;
                phi(l,16) = (y(l,1) - p(l,1)) - rho_val;
            }
        }

		Matrix temp_yp = f(x, y); // Matriz auxiliar
		yp = temp_yp;

        for (int l = 1; l <= n_eqn; l++) {
            phi(l,kp1+1) = yp(l,1) - phi(l,2);
            phi(l,kp2+1) = phi(l,kp1+1) - phi(l,kp2+1);
        }
        for (int i = 1; i <= k; i++) {
            for (int l = 1; l <= n_eqn; ++l) {
                phi(l,i+1) += phi(l,kp1+1);
            }
        }
		
			cout << "17 del dei" << endl;
			
			

        // Estimar error en orden k+1
        erkp1 = 0.0;
        if (knew == km1 || k == 12) {
            phase1 = false;
        }

        if (phase1) {
            k = kp1;
            erk = erkp1;
        } else if (knew == km1) {
            k = km1;
            erk = erkm1;
        } else if (kp1 <= ns) {
            for (int l = 1; l <= n_eqn; l++) {
                erkp1 += std::pow(phi(l,kp2+1) / wt(l,1), 2);
            }
            erkp1 = absh * gstr(kp1+1,1) * std::sqrt(erkp1);
            if (k > 1) {
                if (erkm1 <= std::min(erk, erkp1)) {
                    k = km1;
                    erk = erkm1;
                } else if (erkp1 < erk && k != 12) {
                    k = kp1;
                    erk = erkp1;
                }
            } else if (erkp1 < 0.5 * erk) {
                k = kp1;
                erk = erkp1;
            }
        }


		cout << "18 del dei" << endl;
	
	
	
        // Determinar nuevo tamaño de paso
        if (phase1 || p5eps >= erk * two(k+2,1)) {
				cout << "entras en 18.1 del dei" << endl;
            hnew = 2.0 * h;
        } else if (p5eps < erk) {
			cout << "entras en 18.2 del dei" << endl;
            double temp2 = k + 1;
            double r = std::pow(p5eps / erk, 1.0 / temp2);
            hnew = absh * std::max(0.5, std::min(0.9, r));
            hnew = sign_(std::max(hnew, fouru * std::abs(x)), h);
        } else {
							cout << "entras en 18.3 del dei" << endl;
            hnew = h;
        }
						cout << "entras en 18.4 del dei" << endl;
        h = hnew;

        // Verificar rigidez
        nostep++;
        kle4++;
		
        if (kold > 4) {
							cout << "entras en 18.5 del dei" << endl;
            kle4 = 0;
        }
        if (kle4 >= 50) {
							cout << "entras en 18.6 del dei" << endl;
            stiff = true;
            std::cerr << "PROBELMA DE kle4, kle4=" << kle4 << std::endl;
            State_ = DE_STATE.DE_STIFF;
            y = yy;
            t = x;
            told = t;
            OldPermit = true;
            return y;
        }
    }
	
	cout << "19 del dei, justo antes de return y" << endl;

    return y;
}