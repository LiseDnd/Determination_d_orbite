#include <math.h> 
#include <stdio.h>
#include <string.h>

double a(double x[], double y[], int n); // Fonction coefficient directeur regression linéaire
double b(double x[], double y[], double a, int n); // Fonction ordonnée à l'origine regression linéaire
void cross(double u[3], double v[3], double vect[3]); // Produit vectoriel
double dot(double u[3], double v[3]); // Produit scalaire
double norm(double u[3]); // Norme du vecteur u
double det(double u[3], double v[3], double w[3]); // Déterminant d'une matrice 3x3 où u, v et w sont les vecteurs colonnes
double f(double r, double alpha, double beta, double gamma); // Polynome de degrés 8 pour la méthode de Gauss
double f_prime(double r, double alpha, double beta, double gamma); // Dérivée du polynome 
double Newton(double alpha, double beta, double gamma, double r0, double Delta, double NLoop); // Méthode de Newton pour trouver les racines d'un polynome
void Time_RADEC_Sun(double date, double RADEC[2]); // Calcul de l'ascension droite et declinaison pour le soleil
void velocity(double u[3], double v[3], double w[3], double t1, double t2, double t3, double mu, double v_[3]); // Calcul du vecteur vitesse en connaissant trois points
void orbital_elements(double r[3], double v[3], double mu, double elements[8]); // Calcul des éléments orbitaux à partir de la position et de la vitesse


int main() {
	/* ---------------------------- Coordoonnées Astéroide ----------------------------- */ 
	double ARA1 = 43.537, ADEC1 = -8.7833, ARA2 = 54.420, ADEC2 = -12.074, ARA3 = 64.318, ADEC3 = -15.105; // Coordonnées (ascension droite et declinaison) de l'astéroide
	
	
	
	/* ------------------------------- Méthode de Gauss -------------------------------- */
	double t1 = 0, t2 = 118.10/3600, t3 = 237.58/3600; // Heure d'observation
	double Obs_Lat = +43 + 26/60 +39/3600, Obs_Long = -(01 + 48./60 +58/3600); // Latitude et longitude (en degrés de l'observatoire)
	double HSP = 7.292115e-5; // Mouvement moyen de la Terre rad s-1
	double Rearth = 6378.137e3; // Rayon de la Terre en m
	double Dsun_earth = 149597870.700e3; // Distance Terre/Soleil en m
	double mu_Earth = 398600.441e9/pow(3600, -2); // mu de la Terre m3 h-2 
	
	// Heure sidérale de l'astéroide pour les trois FITS (en radian)
	double HS1 = 44.506*M_PI/180, HS2 = 45.000*M_PI/180, HS3 = 45.499*M_PI/180; 
	
	// Vecteur unitaire observateur/asteroide pour chaque FITS 
	double L1[3] = {cos(ADEC1 *M_PI/180)*cos(ARA1 *M_PI/180), cos(ADEC1 *M_PI/180)*sin(ARA1 *M_PI/180), sin(ADEC1 *M_PI/180)};
	double L2[3] = {cos(ADEC2 *M_PI/180)*cos(ARA2 *M_PI/180), cos(ADEC2 *M_PI/180)*sin(ARA2 *M_PI/180), sin(ADEC2 *M_PI/180)};
	double L3[3] = {cos(ADEC3 *M_PI/180)*cos(ARA3 *M_PI/180), cos(ADEC3 *M_PI/180)*sin(ARA3 *M_PI/180), sin(ADEC3 *M_PI/180)};
	  
	// Vecteur observateur/centre de la Terre pour chaque FITS    
	double R1[3] = {Rearth*cos(Obs_Lat *M_PI/180)*cos(HS1), Rearth*cos(Obs_Lat *M_PI/180)*sin(HS1), Rearth*sin(Obs_Lat *M_PI/180)};                
	double R2[3] = {Rearth*cos(Obs_Lat *M_PI/180)*cos(HS2), Rearth*cos(Obs_Lat *M_PI/180)*sin(HS2), Rearth*sin(Obs_Lat *M_PI/180)};
	double R3[3] = {Rearth*cos(Obs_Lat *M_PI/180)*cos(HS3), Rearth*cos(Obs_Lat *M_PI/180)*sin(HS3), Rearth*sin(Obs_Lat *M_PI/180)};
	
	double Dt1 = t2-t1, Dt2 = t3-t2, Dt3 = t3-t1;
	
	double A1 = Dt1/Dt3, A2 = Dt2/Dt3; 
	double B1 = (mu_Earth*(pow(Dt3, 2)-pow(Dt1, 2))*Dt1)/(6*Dt3), B2 = (mu_Earth*(pow(Dt3, 2)-pow(Dt2, 2))*Dt2)/(6*Dt3);
	
	for (int i = 0; i < 3; i++) {
		// Calcul l'opposé des vecteurs L2 et R2
    		L2[i] = -L2[i];
    		R2[i] = -R2[i];
	}

	double D = det(L1, L2, L3), D1 = det(L1, R1, L3), D2 = det(L1, R2, L3), D3 = det(L1, R3, L3); // Calcul de déterminant de matrices 3x3
	
	double A = D2/D + A2*D1/D + A1*D3/D, B = B2*D1/D + B1*D3/D, C = 2*(R2[0]*L2[0] + R2[1]*L2[1] + R2[2]*L2[2]);
	double alpha = A*C - pow(A,2) - pow(norm(R2), 2), beta = B*C - 2*A*B, gamma = -pow(B,2); // Coefficients polynome en r2
	
	double r2 = Newton(alpha, beta, gamma, 3*Dsun_earth, 1e-9, 1e6); // Norme du vecteur géocentrique de l'astéroide sur le FITS 2 en m
	
	double S13 = (Dt1/Dt3)*(1 + mu_Earth*(pow(Dt3,2)-pow(Dt1,2))/(6*pow(r2,3))), S23 = (Dt2/Dt3)*(1 + mu_Earth*(pow(Dt3,2)-pow(Dt2,2))/(6*pow(r2,3))), S12 = S13/S23; // Rapport des aires
	
	double c1[3], c2[3], c3[3];
	cross(L2, L3, c1);
	cross(L3, L1, c2);
	cross(L1, L2, c3);
	double c123 = dot(c3, L3);
	double c = dot(L1, c1);
	
	for (int i = 0; i < 3; i++) {
    		c1[i] = c1[i]/c123;
    		c2[i] = c2[i]/c123;
    		c3[i] = c3[i]/c123;
	}
	
	for (int i = 0; i < 3; i++) {
		// Calcul l'opposé des vecteurs L2 et R2
    		L2[i] = -L2[i];
    		R2[i] = -R2[i];
	}
	
	double d1 = (-dot(c1, R1) + 1/S23*dot(c1, R2) - S12*dot(c1, R3)), d2 = (-S23*dot(c2, R1) + dot(c2, R2) - S13*dot(c2, R3)), d3 = (-(1/S12)*dot(c3, R1) + 1/S13*dot(c3, R2) - dot(c3, R3));
	
	double vectT_r1[3] = {d1*L1[0]+R1[0], d1*L1[1]+R1[1], d1*L1[2]+R1[2]};
	double vectT_r2[3] = {d2*L2[0]+R2[0], d2*L2[1]+R2[1], d2*L2[2]+R2[2]};
	double vectT_r3[3] = {d3*L3[0]+R3[0], d3*L3[1]+R3[1], d3*L3[2]+R3[2]}; // Vecteur position de l'astéroïde pour les 3 Fits en m (géocentrique)
	
	
	
	/* --------------------- Conversion (r1, r2, r3) -> (r2, v2) ----------------------- */
	double vect_v[3]; // Vecteur vitesse de l'astéroïde
	velocity(vectT_r1, vectT_r2, vectT_r3, t1, t2, t3, mu_Earth, vect_v); // Calcul du vecteur vitesse de l'astéroïde en m h-1
	
	
	
	/* -------------- Conversion (r2, v2) -> (a, e, i, omega, Omega, M) ---------------- */
	double Orbital_Elements[8];
	orbital_elements(vectT_r2, vect_v, mu_Earth, Orbital_Elements);



	/* ---------------------------------- Affichage ------------------------------------ */
	printf("---------------------------- Coordoonnées Astéroide -----------------------------\n");
	printf("FITS1: Ascension Droite: %.9f, Déclinaison: %.9f\n", ARA1, ADEC1);
	printf("FITS2: Ascension Droite: %.9f, Déclinaison: %.9f\n", ARA2, ADEC2);
	printf("FITS3: Ascension Droite: %.9f, Déclinaison: %.9f\n", ARA3, ADEC3);
	printf("\n\n\n");
	
	printf("------------------------------- Méthode de Gauss --------------------------------\n");
	printf("------------------------------- Heure Sidérale \n");
	printf("Heure Sidérale:  Fits 1: %f degrés\n 	         Fits 2: %f degrés\n 	         Fits 3: %f degrés\n", fmod(HS1 *180/M_PI, 360), fmod(HS2 *180/M_PI, 360), fmod(HS3 *180/M_PI, 360));
	printf("\n");
	
	printf("------------------------------- Vecteurs Direction \n");
	printf("FITS1:  L1:");
	for (int i = 0; i < 3; i++) {
    		printf("%f, ", L1[i]);
    		}
    	printf("\n	R1:");
	for (int i = 0; i < 3; i++) {
    		printf("%f, ", R1[i]);
    		}
    	printf("\n");
	printf("FITS2:  L2:");
	for (int i = 0; i < 3; i++) {
    		printf("%f, ", L2[i]);
    		}
    	printf("\n	R2:");
	for (int i = 0; i < 3; i++) {
    		printf("%f, ", R2[i]);
    		}
    	printf("\n");
    	printf("FITS3:  L3:");
	for (int i = 0; i < 3; i++) {
    		printf("%f, ", L3[i]);
    		}
    	printf("\n	R3:");
	for (int i = 0; i < 3; i++) {
    		printf("%f, ", R3[i]);
    		}
    	printf("\n\nNorme de L: Fits 1:%f km\n            Fits 2:%f km\n            Fits 3:%f km", norm(L1)*1e-3, norm(L2)*1e-3, norm(L3)*1e-3);
    	printf("\nNorme de R: Fits 1:%f UA\n            Fits 2:%f UA\n            Fits 3:%f UA", norm(R1)/Dsun_earth, norm(R2)/Dsun_earth, norm(R3)/Dsun_earth);
	printf("\n\n");
	
	printf("------------------------------- Ecarts \n");
	printf("Ecart:  Dt1: %f\n 	Dt2: %f\n 	Dt3: %f\n", Dt1, Dt2, Dt3);
	printf("\n");
	
	printf("------------------------------- Vecteurs géocentrique de l'astéroïde \n");
	printf("Vecteur: r1:");
	for (int i = 0; i < 3; i++) {
    		printf("%f, ", vectT_r1[i]*1e-3);
    		}
	printf("  km   -> r1 = %f km\n", norm(vectT_r1)*1e-3);
	printf("  	 r2:");
	for (int i = 0; i < 3; i++) {
    		printf("%f, ", vectT_r2[i]*1e-3);
    		}
	printf("km   -> r2 = %f km\n", norm(vectT_r2)*1e-3);
	printf("  	 r3:");
	for (int i = 0; i < 3; i++) {
    		printf("%f, ", vectT_r3[i]*1e-3);
    		}
	printf("  km   -> r3 = %f km\n", norm(vectT_r3)*1e-3);
	printf("\n\n\n");
	
	printf("--------------------- Conversion (r1, r2, r3) -> (r2, v2) -----------------------\n");
	printf("------------------------------- Vitesse \n");
	printf("v:");
	for (int i = 0; i < 3; i++) {
    		printf("%f, ", vect_v[i]*1e-3/3600);
    		}
	printf("km s-1\n");
	printf("v=%f km s-1\n", norm(vect_v)*1e-3/3600);
	printf("\n\n\n");
	
	printf("-------------- Conversion (r2, v2) -> (a, e, i, omega, Omega, M) ----------------\n");
	printf("Elements orbitaux: a:     %.9f km \n", Orbital_Elements[0]*1e-3);
	printf("		   e:     %.9f\n", Orbital_Elements[1]);
	printf("		   i:     %f degrés\n", Orbital_Elements[2]);
	printf("		   Omega: %f degrés\n", Orbital_Elements[3]);
	printf("		   omega: %f degrés\n", Orbital_Elements[4]);
	printf("		   Theta: %f degrés\n", Orbital_Elements[5]);
	printf("		   M:     %f degrés\n", Orbital_Elements[6]);
	printf("\n");
}


double a(double x[], double y[], int n) {
	/* ----------------- Fonction coefficient directeur regression linéaire ---------------- */
    	double sum_xy = 0.0, sum_x = 0.0, sum_y = 0.0, sum_x2 = 0.0;

    	for (int i = 0; i < n; i++) {
        	sum_xy += x[i] * y[i];
        	sum_x += x[i];
        	sum_y += y[i];
        	sum_x2 += x[i] * x[i];
    	}

    	double a_num = (n * sum_xy) - (sum_x * sum_y);
    	double a_deno = (n * sum_x2) - (sum_x * sum_x);

    	return a_num / a_deno;
}


double b(double x[], double y[], double a, int n) {
	/* ----------------- Fonction ordonnée à l'origine regression linéaire ----------------- */
    	double sum_x = 0.0, sum_y = 0.0;

    	for (int i = 0; i < n; i++) {
        	sum_x += x[i];
        	sum_y += y[i];
    	}

    	return (sum_y - a * sum_x) / n;
}


void cross(double u[3], double v[3], double vect[3]) {
	/* -------------------------------- Produit vectoriel ---------------------------------- */
	vect[0] = u[1]*v[2] - u[2]*v[1];
	vect[1] = u[2]*v[0] - u[0]*v[2];
	vect[2] = u[0]*v[1] - u[1]*v[0];
}


double dot(double u[3], double v[3]) {
	/* --------------------------------- Produit scalaire ---------------------------------- */
	return u[0]*v[0] + u[1]*v[1] + u[2]*v[2];
}


double norm(double u[3]) {
	/* -------------------------------- Norme d'un vecteur --------------------------------- */
	return sqrt(pow(u[0], 2) + pow(u[1], 2) + pow(u[2], 2));
}


double det(double u[3], double v[3], double w[3]) {
	/* ------- Déterminant d'une matrice 3x3 où u, v et w sont les vecteurs colonnes ------- */
	return u[0]*(v[1]*w[2]-v[2]*w[1]) - u[1]*(v[0]*w[2]-v[2]*w[0]) + u[2]*(v[0]*w[1]-v[1]*w[0]);
}


double f(double r, double alpha, double beta, double gamma) {
	/* ------------------- Polynome de degrés 8 pour la méthode de Gauss ------------------- */
    	return pow(r, 8) + alpha * pow(r, 6) + beta * pow(r, 3) + gamma;
}


double f_prime(double r, double alpha, double beta, double gamma) {
	/* ------------------------------- Dérivée du polynome --------------------------------- */
    	return 8 * pow(r, 7) + 6 * alpha * pow(r, 5) + 3 * beta * pow(r, 2);
}


double Newton(double alpha, double beta, double gamma, double r0, double Delta, double NLoop) {
	/* ------------- Méthode de Newton pour trouver les racines d'un polynome -------------- */
    	double r = r0;
    	double delta;
    	int count = 0;

    	do {
        	delta = f(r, alpha, beta, gamma) / f_prime(r, alpha, beta, gamma);
        	r -= delta;
        	count ++;
        	
        	if (count == NLoop) {
        		break;
        	}
        	
    	} while (fabs(delta) > Delta);

    	return r;
}


void Time_RADEC_Sun(double date, double RADEC[2]) {
	/* ----------- Calcul de l'ascension droite et declinaison pour le soleil -------------- */
	double oblecl_Sun = 23.4393 - 3.563E-7 * date; // Obliquity of the ecliptic en degrés
	double omega_Sun = fmod(282.9404 + 4.70935E-5 * date, 360); // Longitude of perihelion en degrés 
	double a_Sun = 149597870.700e3; // Mean distance en m
	double e_Sun = 0.016709 - 1.151E-9 * date; // Eccentricity
	double M_Sun = fmod(356.0470 + 0.9856002585 * date, 360); // Mean anomaly en degrés
	double L_Sun = fmod(M_Sun + omega_Sun, 360); // Mean londitude en degrés
	double E_Sun = fmod(M_Sun + (180/M_PI) * e_Sun * sin(M_Sun * M_PI/180) * (1.0 + e_Sun * cos(M_Sun *M_PI/180)), 360); // Mean eccentricity en degrés
	double vect_r_Sun[3] = {cos(E_Sun *M_PI/180) - e_Sun, sin(E_Sun * M_PI/180) * sqrt(1 - e_Sun*e_Sun), 0}; // Vecteur coordonnées rectangulaire du soleil dans le plan de l'ecliptique (x pointant vers la périhélie)
	double r_Sun = norm(vect_r_Sun); // Distance Terre/Soleil dans le plan de l'ecliptique
    	double nu_Sun = atan2(vect_r_Sun[1], vect_r_Sun[0])*180/M_PI; // Anomalie vraie en degrés
    	double Long_Sun = fmod(nu_Sun + omega_Sun, 360); // Longitude du soleil en degrés
    	double r_ecli[3] = {r_Sun * cos(Long_Sun *M_PI/180), r_Sun * sin(Long_Sun *M_PI/180), 0.0}; // Vecteur coordonnées du plan de l'ecliptique 
    	double r_equa[3] = {r_ecli[0], r_ecli[1] * cos(oblecl_Sun *M_PI/180), r_ecli[1] * sin(oblecl_Sun *M_PI/180)}; // Vecteur coordonnées du plan equatorial ???????????????????????????????????????????????
    	double r_Sun_bis = norm(r_equa); // Distance Terre/Soleil dans le plan equatorial
    	double Sun_RA =  fmod(atan2(r_equa[1], r_equa[0]) *180/M_PI, 360); // Ascension droite du Soleil
    	double Sun_DEC =  fmod(atan2(r_equa[2], sqrt(r_equa[0]*r_equa[0] + r_equa[1]*r_equa[1])) *180/M_PI, 360); // Declinaison du Soleil
    	
    	RADEC[0] = Sun_RA;
    	RADEC[1] = Sun_DEC;
}


void velocity(double u[3], double v[3], double w[3], double t1, double t2, double t3, double mu, double v_[3]) {
	/* -------------- Calcul du vecteur vitesse en connaissant trois points ---------------- */
	double tau1 = t1 - t2, tau3 = t3-t2;
	double f1 = 1 - (1/2)*(mu/pow(norm(v), 3))*pow(tau1, 2), f3 = 1 - (1/2)*(mu/pow(norm(v), 3))*pow(tau3, 2);
	double g1 = tau1 - (1/6)*(mu/pow(norm(v), 3))*pow(tau1, 3), g3 = tau3 - (1/6)*(mu/pow(norm(v), 3))*pow(tau3, 3);
	
	for (int i = 0; i < 3; i++) {
		v_[i] = (-f3*u[i] + f1*w[i])/(f1*g3 - f3*g1);
	}
}


void orbital_elements(double r[3], double v[3], double mu, double elements[8]) {
	/* ------ Calcul des éléments orbitaux à partir de la position et de la vitesse -------- */
	double r_norm = norm(r), v_norm = norm(v);
	double vr = dot(r, v)/r_norm;
	double h[3], n[3];
	double Omega, omega;
	double Z[3] = {0, 0, 1};
	
	cross(r, v, h);
	cross(Z, h, n);

	double h_norm = norm(h), n_norm = norm(n);

	double i = acos(h[2]/h_norm);
	
	Omega = acos(n[0]/n_norm);
	if (n[1] < 0) {
		Omega = 2*M_PI - Omega;
	}
		
	double e[3] = {1/mu*((pow(v_norm, 2) - mu/r_norm)*r[0] - r_norm*vr*v[0]), 1/mu*((pow(v_norm, 2) - mu/r_norm)*r[1] - r_norm*vr*v[1]), 1/mu*((pow(v_norm, 2) - mu/r_norm)*r[2] - r_norm*vr*v[2])};
	double e_norm = norm(e);
	
	omega = acos(dot(n, e)/(n_norm*e_norm));
	
	if (e[2] < 0) {
		omega = 2*M_PI -omega;
	}

	double TA = acos(dot(e,r)/(e_norm*r_norm));
	if (vr < 0) {
		TA = 2*M_PI - TA;
	}
	
	double E = 2*atan2(tan(TA/2), sqrt((1+e_norm)/(1-e_norm)));
	double M = E - e_norm*sin(E);
	
	double a = pow(h_norm, 2)/(mu*(1 - pow(e_norm, 2)));
	
	elements[0] = a;
	elements[1] = e_norm;
	elements[2] = i *180/M_PI;
	elements[3] = Omega *180/M_PI;
	elements[4] = omega *180/M_PI;
	elements[5] = TA *180/M_PI;
	elements[6] = M *180/M_PI;
	elements[7] = h_norm *180/M_PI;
}

