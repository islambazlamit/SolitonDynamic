# SolitonDynamic
This code is designed to study soliton dynamics by solving the 1D nonlinear Schrödinger equation (NLSE) with dual nonlinearity and an external potential, using the Iterative Power Series (IPS) method.
The original implementation, written in FORTRAN by Prof Usama AL Khawaja, provided 16-digit precision, whereas this C++ version extends the precision significantly, offering up to 50 to 100 digits of accuracy. 

                                          iψt​(x,t)+a1​ψxx​(x,t)+a2​∣ψ(x,t)∣2ψ(x,t)+a3​∣ψ(x,t)∣4ψ(x,t)+V(x)ψ(x,t)=0

You may need to choose your solution and potential. It is valid for single and two-soliton solutions. 


This version is created by Islam M. Albazlamit, islam.albazlamit@outlook.com
