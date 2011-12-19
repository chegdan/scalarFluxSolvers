//determine the initial values of Dturbulent

if(tSchmidt>0){//if this is turned on in the SIMPLE area of fvSolutions
//read and determine the Dturbulent from the turbulent schmidt number and turbulent viscosity
	Info<<"\nDetermining Dturbulence from "<<turbulentSchmidt<<endl;
	Dturbulent = nut/turbulentSchmidt;
}else {
	Info<<"\nDturbulence = 0 over domain"<<endl;
	Dturbulent = D-D;//essentially zero
	}
