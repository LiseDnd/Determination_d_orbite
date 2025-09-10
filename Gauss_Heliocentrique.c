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
void orbital_elements(double r[3], double v[3], double date, double mu, double elements[8]); // Calcul des éléments orbitaux à partir de la position et de la vitesse


int main() {
	/* ----------------------------- Regression Linéaire ------------------------------- */ 
	int n = 10; // nombre de points pour la regression linéaire
	double SRA[10] = {84.973940, 84.957663, 84.952060, 84.939592, 85.012794, 84.937589, 84.922586, 84.968012, 84.956148, 84.954905}, 
	      SDEC[10] = {28.484180, 28.489463, 28.492225, 28.498026, 28.450026, 28.484975, 28.461919, 28.452122, 28.442803, 28.445037}; 
	      // Coordonnées en ascension droite et déclinaison de 10 étoiles repérées sur les images (en dregrés)
	double SPA1[10] = {2141.47, 2228.01, 2257.87, 2324.43, 1933.45, 2334.40, 2411.88, 2170.61, 2232.88, 2239.49}, SPO1[10] = {1311.36, 1280.91, 1264.92, 1230.51, 1514.42, 1309.39, 1449.23, 1504.77, 1561.91, 1548.58}, 
	      SPA2[10] = {1877.7, 1964.2, 1994.1, 2061.0, 1669.1, 2070.2, 2147.7, 1907.7, 1969.0, 1976.6}, SPO2[10] = {1489.5, 1458.3, 1442.2, 1408.8, 1692.3, 1487.7, 1627.1, 1683.2, 1740.3, 1726.6}, 
	      SPA3[10] = {1766.6, 1854.53, 1884.58, 1950.37, 1559.73, 1960.27, 2037.95, 1796.65, 1858.68, 1865.60}, SPO3[10] = {1528.7, 1498.52, 1482.74, 1448.01, 1731.85, 1527.06, 1666.80, 1722.51, 1779.35, 1766.49}; 
	      // Coordonnées en pixels de 10 étoiles repérées sur les images pour les trois images
	
	double a1_RA = a(SPA1, SRA, n), a1_DEC = a(SPO1, SDEC, n), a2_RA = a(SPA2, SRA, n), a2_DEC = a(SPO2, SDEC, n), a3_RA = a(SPA3, SRA, n), a3_DEC = a(SPO3, SDEC, n); // Coefficient directeur regression linéaire pour les trois FITS (abcisse et ordonnée)
	double b1_RA = b(SPA1, SRA, a1_RA, n), b1_DEC = b(SPO1, SDEC, a1_DEC, n), b2_RA = b(SPA2, SRA, a2_RA, n), b2_DEC = b(SPO2, SDEC, a2_DEC, n), b3_RA = b(SPA3, SRA, a3_RA, n), b3_DEC = b(SPO3, SDEC, a3_DEC, n); // Ordonnée à l'origine regression linéaire pour les trois FITS (abcisse et ordonnée)
	
	
	 
	/* ---------------------------- Coordoonnées Astéroide ----------------------------- */ 
	double APA1 = 2290.05, APO1 = 1407.82, APA2 = 1954.16, APO2 = 1577.93, APA3 = 1776.49, APO3 = 1609.38; // Coordonnées (pixels) de l'astéroide sur les trois FITS
	double ARA1 = a1_RA*APA1 + b1_RA, ADEC1 = a1_DEC*APO1 + b1_DEC, ARA2 = a2_RA*APA2 + b2_RA, ADEC2 = a2_DEC*APO2 + b2_DEC, ARA3 = a3_RA*APA3 + b3_RA, ADEC3 = a3_DEC*APO3 + b3_DEC; 
	      // Coordonnées (ascension droite et declinaison) de l'astéroide sur les trois FITS
	
	
	
	/* ----------------------------- Coordoonnées Soleil ------------------------------- */ 
	double t1 = 19 + 40./60 + 3.44/3600, t2 = 20 + 59./60 + 34.41/3600, t3 = 22 +12./60 + 37.57/3600; // Heure d'observation pour les trois FITS
	double DN1_Sun = 367*2008 - 7 * ( 2008 + (3+9)/12 ) / 4 - 3 * ( ( 2008 + (3-9)/7 ) / 100 + 1 ) / 4 + 275*3/9 + 2 - 730515 + t1/24, DN2_Sun = 367*2008 - 7 * ( 2008 + (3+9)/12 ) / 4 - 3 * ( ( 2008 + (3-9)/7 ) / 100 + 1 ) / 4 + 275*3/9 + 2 - 730515 + t2/24, DN3_Sun = 367*2008 - 7 * ( 2008 + (3+9)/12 ) / 4 - 3 * ( ( 2008 + (3-9)/7 ) / 100 + 1 ) / 4 + 275*3/9 + 2 - 730515 + t3/24; 
	
	double Sun_RADEC1[2], Sun_RADEC2[2], Sun_RADEC3[2];
	
	Time_RADEC_Sun(DN1_Sun, Sun_RADEC1);
	Time_RADEC_Sun(DN2_Sun, Sun_RADEC2);
	Time_RADEC_Sun(DN3_Sun, Sun_RADEC3);
	
	
	
	/* ------------------------------- Méthode de Gauss -------------------------------- */
	double Obs_Lat = +43 + 26/60 +39/3600, Obs_Long = -(01 + 48./60 +58/3600); // Latitude et longitude (en degrés de l'observatoire)
	double HSP = 7.292115e-5; // Mouvement moyen de la Terre rad s-1
	double Rearth = 6378.137e3; // Rayon de la Terre en m
	double Dsun_earth = 149597870.700e3; // Distance Terre/Soleil en m
	double mu_Sun = 132712440018e9/pow(3600, -2), mu_Earth = 398600.441e9/pow(3600, -2); // mu du Soleil m3 h-2 et mu de la Terre m3 h-2 
	
	// Vecteur unitaire observateur/asteroide pour chaque FITS 
	double L1[3] = {cos(ADEC1 *M_PI/180)*cos(ARA1 *M_PI/180), cos(ADEC1 *M_PI/180)*sin(ARA1 *M_PI/180), sin(ADEC1 *M_PI/180)};
	double L2[3] = {cos(ADEC2 *M_PI/180)*cos(ARA2 *M_PI/180), cos(ADEC2 *M_PI/180)*sin(ARA2 *M_PI/180), sin(ADEC2 *M_PI/180)};
	double L3[3] = {cos(ADEC3 *M_PI/180)*cos(ARA3 *M_PI/180), cos(ADEC3 *M_PI/180)*sin(ARA3 *M_PI/180), sin(ADEC3 *M_PI/180)};
	
	// Vecteur unitaire Soleil/Terre
	double S1[3] = {-Dsun_earth*cos(Sun_RADEC1[1] *M_PI/180)*cos(Sun_RADEC1[0] *M_PI/180), -Dsun_earth*cos(Sun_RADEC1[1] *M_PI/180)*sin(Sun_RADEC1[0] *M_PI/180), -Dsun_earth*sin(Sun_RADEC1[1] *M_PI/180)},
	       S2[3] = {-Dsun_earth*cos(Sun_RADEC2[1] *M_PI/180)*cos(Sun_RADEC2[0] *M_PI/180), -Dsun_earth*cos(Sun_RADEC2[1] *M_PI/180)*sin(Sun_RADEC2[0] *M_PI/180), -Dsun_earth*sin(Sun_RADEC2[1] *M_PI/180)},
	       S3[3] = {-Dsun_earth*cos(Sun_RADEC3[1] *M_PI/180)*cos(Sun_RADEC3[0] *M_PI/180), -Dsun_earth*cos(Sun_RADEC3[1] *M_PI/180)*sin(Sun_RADEC3[0] *M_PI/180), -Dsun_earth*sin(Sun_RADEC3[1] *M_PI/180)}; 
	
	double Dt1 = t2-t1, Dt2 = t3-t2, Dt3 = t3-t1;
	
	double A1 = Dt1/Dt3, A2 = Dt2/Dt3; 
	double B1 = (mu_Sun*(pow(Dt3, 2)-pow(Dt1, 2))*Dt1)/(6*Dt3), B2 = (mu_Sun*(pow(Dt3, 2)-pow(Dt2, 2))*Dt2)/(6*Dt3);
	
	for (int i = 0; i < 3; i++) {
		// Calcul l'opposé des vecteurs L2 et R2
    		L2[i] = -L2[i];
    		S2[i] = -S2[i];
	}

	double D = det(L1, L2, L3), D1 = det(L1, S1, L3), D2 = det(L1, S2, L3), D3 = det(L1, S3, L3); // Calcul de déterminant de matrices 3x3
	
	double A = D2/D + A2*D1/D + A1*D3/D, B = B2*D1/D + B1*D3/D, C = 2*(S2[0]*L2[0] + S2[1]*L2[1] + S2[2]*L2[2]);
	double alpha = A*C - pow(A,2) - pow(norm(S2), 2), beta = B*C - 2*A*B, gamma = -pow(B,2); // Coefficients polynome en r2
	
	double r2 = Newton(alpha, beta, gamma, 3*Dsun_earth, 1e-9, 1e6); // Norme du vecteur géocentrique de l'astéroide sur le FITS 2 en m
	
	double S13 = (Dt1/Dt3)*(1 + mu_Sun*(pow(Dt3,2)-pow(Dt1,2))/(6*pow(r2,3))), S23 = (Dt2/Dt3)*(1 + mu_Sun*(pow(Dt3,2)-pow(Dt2,2))/(6*pow(r2,3))), S12 = S13/S23; // Rapport des aires
	
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
    		S2[i] = -S2[i];
	}
	
	double d1 = (-dot(c1, S1) + 1/S23*dot(c1, S2) - S12*dot(c1, S3)), d2 = (-S23*dot(c2, S1) + dot(c2, S2) - S13*dot(c2, S3)), d3 = (-(1/S12)*dot(c3, S1) + 1/S13*dot(c3, S2) - dot(c3, S3));
	
	double vect_r1[3] = {d1*L1[0]+S1[0], d1*L1[1]+S1[1], d1*L1[2]+S1[2]};
	double vect_r2[3] = {d2*L2[0]+S2[0], d2*L2[1]+S2[1], d2*L2[2]+S2[2]};
	double vect_r3[3] = {d3*L3[0]+S3[0], d3*L3[1]+S3[1], d3*L3[2]+S3[2]}; // Vecteur position de l'astéroïde pour les 3 Fits en m (héliocentrique)
	
	
	
	/* --------------------- Conversion (r1, r2, r3) -> (r2, v2) ----------------------- */
	double vect_v[3]; // Vecteur vitesse de l'astéroïde
	velocity(vect_r1, vect_r2, vect_r3, t1, t2, t3, mu_Sun, vect_v); // Calcul du vecteur vitesse de l'astéroïde en m h-1
	
	
	
	/* -------------- Conversion (r2, v2) -> (a, e, i, omega, Omega, M) ---------------- */
	double Orbital_Elements[8];
	orbital_elements(vect_r2, vect_v, DN2_Sun, mu_Sun, Orbital_Elements);
	
	
	
	
	printf("----------------------------- Regression Linéaire -------------------------------\n");
	printf("FITS1: a1_RA: %.9f, a1_DEC: %.9f, b1_RA: %f, b1_DEC: %f\n", a1_RA, a1_DEC, b1_RA, b1_DEC);
	printf("FITS2: a2_RA: %.9f, a2_DEC: %.9f, b2_RA: %f, b2_DEC: %f\n", a2_RA, a2_DEC, b2_RA, b2_DEC);
	printf("FITS3: a3_RA: %.9f, a3_DEC: %.9f, b3_RA: %f, b3_DEC: %f\n", a3_RA, a3_DEC, b3_RA, b3_DEC); 
	printf("\n\n\n");
	
	printf("---------------------------- Coordoonnées Astéroide -----------------------------\n");
	printf("FITS1: Ascension Droite: %.9f, Déclinaison: %.9f\n", ARA1, ADEC1);
	printf("FITS2: Ascension Droite: %.9f, Déclinaison: %.9f\n", ARA2, ADEC2);
	printf("FITS3: Ascension Droite: %.9f, Déclinaison: %.9f\n", ARA3, ADEC3);
	printf("\n\n\n");
	
	printf("----------------------------- Coordoonnées Soleil -------------------------------\n");
	printf("FITS1: Ascension droite: %.9f\n", Sun_RADEC1[0]);
	printf("       Declinaison:      %.9f\n", Sun_RADEC1[1]);
	printf("FITS2: Ascension droite: %.9f\n", Sun_RADEC2[0]);
	printf("       Declinaison:      %.9f\n", Sun_RADEC2[1]);
	printf("FITS3: Ascension droite: %.9f\n", Sun_RADEC3[0]);
	printf("       Declinaison:      %.9f\n", Sun_RADEC3[1]);
	printf("\n\n\n");
	
	printf("------------------------------- Méthode de Gauss --------------------------------\n");
	printf("------------------------------- Vecteurs Direction \n");
	printf("FITS1:  L1:");
	for (int i = 0; i < 3; i++) {
    		printf("%f, ", L1[i]);
    		}
        printf("\n	S1:");
	for (int i = 0; i < 3; i++) {
    		printf("%.9e, ", S1[i]);
    		}
    	printf("\n");
	printf("FITS2:  L2:");
	for (int i = 0; i < 3; i++) {
    		printf("%f, ", L2[i]);
    		}
    	printf("\n	S2:");
	for (int i = 0; i < 3; i++) {
    		printf("%.9e, ", S2[i]);
    		}
    	printf("\n");
    	printf("FITS3:  L3:");
	for (int i = 0; i < 3; i++) {
    		printf("%f, ", L3[i]);
    		}
    	printf("\n	S3:");
	for (int i = 0; i < 3; i++) {
    		printf("%.9e, ", S3[i]);
    		}
    	printf("\n\nNorme de L: Fits 1:%f km\n            Fits 2:%f km\n            Fits 3:%f km", norm(L1)*1e-3, norm(L2)*1e-3, norm(L3)*1e-3);
    	printf("\nNorme de S: Fits 1:%f UA\n            Fits 2:%f UA\n            Fits 3:%f UA", norm(S1)/Dsun_earth, norm(S2)/Dsun_earth, norm(S3)/Dsun_earth);
	printf("\n\n");
	
	printf("------------------------------- Ecarts \n");
	printf("Ecart:  Dt1: %f\n 	Dt2: %f\n 	Dt3: %f\n", Dt1, Dt2, Dt3);
	printf("\n");
	
	printf("------------------------------- Vecteurs héliocentrique de l'astéroïde \n");
	printf("Vecteur: r1:");
	for (int i = 0; i < 3; i++) {
    		printf("%f, ", vect_r1[i]*1e-3);
    		}
	printf("  km   -> r1 = %f UA\n", norm(vect_r1)/Dsun_earth);
	printf("  	 r2:");
	for (int i = 0; i < 3; i++) {
    		printf("%f, ", vect_r2[i]*1e-3);
    		}
	printf("km   -> r2 = %f UA\n", norm(vect_r2)/Dsun_earth);
	printf("  	 r3:");
	for (int i = 0; i < 3; i++) {
    		printf("%f, ", vect_r3[i]*1e-3);
    		}
	printf("  km   -> r3 = %f UA\n", norm(vect_r3)/Dsun_earth);
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
	printf("Elements orbitaux: a:     %.9f km = %.9f UA\n", Orbital_Elements[0]*1e-3, Orbital_Elements[0]/Dsun_earth);
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
    	double r_equa[3] = {r_ecli[0], r_ecli[1] * cos(oblecl_Sun *M_PI/180), r_ecli[1] * sin(oblecl_Sun *M_PI/180)}; // Vecteur coordonnées du plan equatorial
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


void orbital_elements(double r[3], double v[3], double date, double mu, double elements[8]) {
	/* ------ Calcul des éléments orbitaux à partir de la position et de la vitesse -------- */
	double oblecl_Sun = 23.4393 - 3.563E-7 * date; // Obliquity of the ecliptic en degrés
	double Px[3] = {1, 0, 0}, Py[3] = {0, cos(oblecl_Sun *M_PI/180), sin(oblecl_Sun *M_PI/180)}, Pz[3] = {0, -sin(oblecl_Sun *M_PI/180), cos(oblecl_Sun *M_PI/180)};
	
	r[0] = dot(r, Px), r[1] = dot(r, Py), r[2] = dot(r, Pz);
	v[0] = dot(v, Px), v[1] = dot(v, Py), v[2] = dot(v, Pz);
	
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

