//determine the initial values of Dturbulent

if(tSchmidt>0){//if this is turned on in the SIMPLE area of fvSolutions
//read and determine the Dturbulent from the turbulent schmidt number and turbulent viscosity
	Info<<"\nDetermining Dturbulence from "<<Sct<<endl;
	Dt = nut/Sct;
}else {
	Info<<"\nDturbulence = 0 over domain"<<endl;
	Dt = D-D;//essentially zero
	}
