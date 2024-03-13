// Created by Alex O'Hare, 13/02/2024
// Hodgkin-Huxley, attempt #1

// Import required packages
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <iomanip>

// Define current as a global variable
float Iext = 0.1 / (7.854 * pow(10, -3));

// Define the Hodkin-Huxley differential equations
void hh(float t, float V, float m, float h, float n, float & alpha_m, float & beta_m, float & alpha_h, float & beta_h, float & alpha_n, float & beta_n, float & dVdt, float & dmdt, float & dhdt, float & dndt){
	const float Cm = 1.0, VNa = 50.0, VK = -77.0, VL = -54.387, gNa = 120.0, gK = 36.0, gL = 0.3;
	alpha_m = (0.1 * (-V - 40.)) / (exp((-V - 40.) / 10.) - 1.);
	beta_m =  4. * exp((-V - 65.) / 18.);
	alpha_h = 0.07 * exp((-V - 65.) / 18.);
	beta_h = 1. / (exp((-V - 35.) / 10.) + 1.);
	alpha_n = (0.01 * (-V - 55.)) / (exp((-V - 55.) / 10.) -1.);
	beta_n = 0.125 * exp((-V - 65.) / 80.); 
	dVdt = (1. / Cm) * (Iext - gL * (V - VL) - gNa * m * m * m * h * (V - VNa) - gK * n * n * n * n * (V - VK));
	dmdt = (alpha_m * (1 - m)) - (beta_m * m);
	dhdt = (alpha_h * (1 - h)) - (beta_h * h);
	dndt = (alpha_n * (1 - n)) - (beta_n * n);}

int main(void){

	// Define initial conditions
	float V = -65.0;
	float m = 0.053;
	float h = 0.6;
	float n = 0.318;

	// Define conditions for Euler method
	float dt = 0.001; // time step

	// Open files to write the data
	std::ofstream outputFile1("hhDataOne.txt");  // V against t
	std::ofstream outputFile2("hhDataTwo.txt"); // m, h, n against t
	
	// Time integration loop using Euler's method
	for(float t = 0; t <= 100; t += dt){
		
		if ( t <= 50) Iext = 0.1 / (7.854 * pow(10, -3)); //0;
		// Output solutions for the file
		outputFile1 << t << "\t" << std::fixed << std::setprecision(4) << V << std::endl;
		outputFile2 << t << "\t" << m << "\t" << h << "\t" << n << std::endl;
		
		if ( t > 50) Iext = 0;
		// Output solutions for the file
		outputFile1 << t << "\t" << std::fixed << std::setprecision(4) << V << std::endl;
		outputFile2 << t << "\t" << m << "\t" << h << "\t" << n << std::endl;
		
		// Calculate the derivatives
		float alpha_m, beta_m, alpha_h, beta_h, alpha_n, beta_n, dVdt, dmdt, dhdt, dndt;
		hh(t, V, m, h, n, alpha_m, beta_m, alpha_h, beta_h, alpha_n, beta_n, dVdt, dmdt, dhdt, dndt);
		
		// Update the solutions using Euler's method
		V += dt * dVdt;
		m += dt * dmdt;
		h += dt * dhdt;
		n += dt * dndt;}
	
	// Close the output files
	outputFile1.close();
	outputFile2.close();
	
	// Use gnuplot to plot the results 
	
	// Plot V vs t
	FILE * gnuplotPipe1 = popen("gnuplot -persistent", "w");
	fprintf(gnuplotPipe1, "plot 'hhDataOne.txt' using 1:2 with lines linecolor 'black' title 'V(t)'\n"); 
	fflush(gnuplotPipe1);
	pclose(gnuplotPipe1);
	
	// Plot m, h, n vs t
	FILE * gnuplotPipe2 = popen("gnuplot -persistent", "w");
	fprintf(gnuplotPipe2, "plot 'hhDataTwo.txt' using 1:2 with lines title 'm(t)', \
								'hhDataTwo.txt' using 1:3 with lines title 'h(t)', \
								'hhDataTwo.txt' using 1:4 with lines title 'n(t)'\n"); 
	fflush(gnuplotPipe2);
	pclose(gnuplotPipe2);
	
	std::cout << "Hodgkin-Huxley simulation completed, data saved in 'hhDataOne.txt' for V vs t and 'hhDataTwo.txt' for m, h, n vs t. Respective plots displayed in gnuplot." << std::endl;
	
	return 0;}
