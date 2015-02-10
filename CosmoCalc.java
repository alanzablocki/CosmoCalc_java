/** CosmoCalc in JAVA V 0.9.4 : Alan Zablocki 2015
 * Code to calculate distances to redshift zmax, in cosmologies where the universe is not flat.
 * Includes radiation and allows dark energy equation of state parameter wde0 differ from -1.
 * Includes PPF parametrization where w(a) = wde0 + wa*(1-a)
 *
 * javac CosmoCalc.java
 * java CosmoCalc
 */

import java.text.DecimalFormat;

public class CosmoCalc
{
	public static DecimalFormat fmt = new DecimalFormat("0.################");
	public static DecimalFormat fmt2 = new DecimalFormat("0.########");

	public static void main(String[] args)
	{

		double tolerance = 1E-6;				// At what point is omk close to == 0 and answers do not change!

		// Define some astronomical constants
		double pc = 3.085677581491164E+16;			// Distance unit, 1 parsec (meters) 
		double SpeedLight = 299792458;				// Speed of light c, in m/s

		/** These quantities will change depending user's choice
		 */

		double littleH = 0.7;
		double HubbleConstant = littleH*100000;			// Hubble expansion rate in m/s/Mpc	
		double omegam = 0.3;					// Matter density parameter
		double omegal = 0.7;					// Dark energy density parameter
		double wde0   = -1.0;					// Equation of state of dark energy
		double wa     = 0.0;					// if w.ne.-1 use PPF with w(a) = wde0 + wa*(1-a)

		double[] params =new double[11];			// Here store cosmological parameters and constants

		params[0] = (2.47240029E-5)/(littleH*littleH);	// Radiation density in photons
		params[1] = omegam;
		params[2] = omegal;					// dark energy
		params[3] = wde0;					// EOS w0
		params[4] = wa;						// EOS wa
		params[5] = HubbleConstant;				// hubble constant 100*h
		params[6] = SpeedLight;					// speed of light
		params[7] = pc*1.0E6;					// Mpc - ALL UNITS ARE CONSISTENT
		params[8] = 1.0 - params[2] - params[1];		// curvature - omk = 1 - omegam - omegal
		params[9] = params[6]/params[5];			// d_H - hubble distance factor
		params[10] = tolerance;
		double factor = (2.0*Math.PI)*Math.pow((params[9]),3.0);

		double zmin  = 0;					// define integration limits, zmin = 0
		double zmax = 0.1;					// zmax is the redshift at which to compute various quantities

final long startTime = System.currentTimeMillis();

		double ComDist = doSimpson(zmin, zmax, params);		// Comoving Distance to zmax
		double age_L   = doSimpsonage(zmin, zmax, params);	// lookback time to zmax
		double age_z0  = doSimpsonage(zmin, 300, params);	// use zmax = 300 to get approximate age of the universe
final long endTime = System.currentTimeMillis();

System.out.println("Calls time: " + (endTime - startTime) );
		double age_Z = age_z0-age_L;				/** Age of the universe at zmax
		 * Total age = lookback time[to zmax] + age[at zmax]
		 */

		//if (Math.abs(params[8]) < tolerance){
		//System.out.println("Universe is FLAT");
		//} else {
		//System.out.println("Universe is not FLAT");
		//}

		//Output results !!!
final long startTime2 = System.currentTimeMillis();

		if (Math.abs(params[8]) < tolerance){

			System.out.println("Universe is FLAT");

			double ComVolflat = (1E-9)*(4.0/3.0)*(Math.PI)*Math.pow(ComDist,3.0); // only true in flat universe

			System.out.println("Comoving Distance Dc at z="+zmax+" is "+
			(double)Math.round(ComDist * 1000) / 1000+" Mpc");
			System.out.println("Angular Diameter Distance DA at z="+zmax+" is "+
			(double)Math.round(ComDist/(1 + zmax) * 1000) / 1000+" Mpc");
			System.out.println("Luminosity Distance DL at z="+zmax+" is "+
			(double)Math.round(ComDist*(1 + zmax) * 1000) / 1000+" Mpc");
			System.out.println("Comoving Volume at z="+zmax+" is "+
			(double)Math.round( ComVolflat * 1000) / 1000+" Gpc^3");
		} else if (params[8] > tolerance){

			System.out.println("Universe is OPEN");

			// get DM
			double coef22 = (Math.pow(Math.abs(params[8]),-0.5))*(params[9]);
			double ans21 = coef22*Math.sinh(Math.pow(Math.abs(params[8]),0.5)*ComDist);
			// get DA and DL
			double ans22 = ans21/(1 + zmax);
			double ans222 = ans21*(1 + zmax);
			// get volume
			double dmbydh = ans21/params[9];
			double dmbydh2 = Math.pow(dmbydh,2);
			double ComVolopen = (1E-9)*(factor/params[8])*
			       (dmbydh*Math.pow((1.0+params[8]*dmbydh2),0.5)-
			       (Math.pow(Math.abs(params[8]),-0.5))*asinh(Math.pow(Math.abs(params[8]),0.5) * dmbydh));

			System.out.println("Comoving Distance Dc at z="+zmax+" is "+
			(double)Math.round(ComDist*(params[9]) * 1000) / 1000+" Mpc");
			System.out.println("Comoving (Transverse) Distance DM at z="+zmax+" is "+
			(double)Math.round(ans21 * 1000) / 1000+" Mpc");
			System.out.println("Angular Diameter Distance DA at z="+zmax+" is "+
			(double)Math.round(ans22 * 1000) / 1000+" Mpc");
			System.out.println("Luminosity Distance DL at z="+zmax+" is "+
			(double)Math.round(ans222 * 1000) / 1000+" Mpc");
			System.out.println("Comoving Volume at z="+zmax+" is "+
			(double)Math.round(ComVolopen * 1000) / 1000+" Gpc^3");

		} else if (params[8] < -tolerance){
			System.out.println("Universe is CLOSED");

			double coef22 = (Math.pow(Math.abs(params[8]),-0.5))*(params[9]);
			double ans21 = coef22*Math.sin(Math.pow(Math.abs(params[8]),0.5)*ComDist);
			double ans22 = ans21/(1 + zmax);
			double ans222 = ans21*(1 + zmax);
			double dmbydh = ans21/params[9];
			double dmbydh2 = Math.pow(dmbydh,2);
			double ComVolclosed = (1E-9)*(factor/params[8])*
			       (dmbydh*Math.pow((1.0+params[8]*dmbydh2),0.5)-
			       (Math.pow(Math.abs(params[8]),-0.5))*Math.asin(Math.pow(Math.abs(params[8]),0.5) * dmbydh));

			System.out.println("Comoving Distance DC at z="+zmax+" is "+
			(double)Math.round(ComDist*(params[9]) * 1000) / 1000+" Mpc");
			System.out.println("Comoving (Transverse) Distance DM at z="+zmax+" is "+
			(double)Math.round(ans21 * 1000) / 1000+" Mpc");
			System.out.println("Angular Diameter Distance DA at z="+zmax+" is "+
			(double)Math.round(ans22 * 1000) / 1000+" Mpc");
			System.out.println("Luminosity Distance DL at z="+zmax+" is "+
			(double)Math.round(ans222 * 1000) / 1000+" Mpc");
			System.out.println("Comoving Volume at z="+zmax+" is "+
			(double)Math.round(ComVolclosed * 1000) / 1000+" Gpc^3");

		} else{
		}

		System.out.println("Age of the Universe at z="+zmax+" is "+
		(double)Math.round(age_Z * 1000) / 1000+" Gyrs");
		System.out.println("Age of the universe today is "+
		(double)Math.round(age_z0 * 1000) / 1000+" Gyrs");
		System.out.println("Lookback time at z="+zmax+" is "+
		(double)Math.round(age_L * 1000) / 1000+" Gyrs");
	
final long endTime2 = System.currentTimeMillis();
System.out.println("output time: " + (endTime2 - startTime2) );

}

	static double asinh(double x){ 
		return Math.log(x + Math.sqrt(x*x + 1.0)); 
	} 

	static double func(double x, double params[])
	{

		if (Math.abs(params[8]) < params[10]){
			return params[9]*Math.pow((params[1]*Math.pow((1.0+x),3.0)+
			params[0]*Math.pow((1.0+x),4.0)+params[2]*Math.pow((1.0+x),3.0*(1.0+params[3]))),-0.5);

		}else{

			return Math.pow((params[1]*Math.pow((1.0+x),3.0)+
			params[8]*Math.pow((1.0+x),2.0)+params[0]*Math.pow((1.0+x),4.0)+
			params[2]*Math.pow((1.0+x),3.0*(1.0+params[3]))),-0.5);
		}
	}

	static double funcage(double x, double params[])
	{       
		double coef3 = (params[7]/params[5])/(31557600.0E9); // get age in Gyrs

		if (Math.abs(params[8]) < params[10]){

			return coef3*(1.0/(1.0+x))*Math.pow((params[1]*Math.pow((1.0+x),3.0)+
			params[0]*Math.pow((1.0+x),4.0)+params[2]*Math.pow((1.0+x),3.0*(1.0+params[3]))),-0.5);
		}else{

			return coef3*(1.0/(1.0+x))*Math.pow((params[1]*Math.pow((1.0+x),3.0)+
			params[8]*Math.pow((1.0+x),2.0)+params[0]*Math.pow((1.0+x),4.0)+
			params[2]*Math.pow((1.0+x),3.0*(1.0+params[3]))),-0.5);
		}
	}

	// ********** Simpson routines **********

	static int    JMAX=20; 
	static double TOL=1E-9;

	static double doSimpson(double a, double b, double fp[])
	{
		double t, tprev=-999;					// trapezoid values
		double p, pprev=-999;					// predicted values
		int nnext = 1;
		t = 0.5*(b-a)*(func(a,fp) + func(b,fp));		// initial guess
		for (int j=1; j<JMAX; j++)
		{
			double delta = (b-a)/nnext;
			double x = a + 0.5*delta;
			double sum = 0.0;					// the trapezoid
			for (int k=0; k<nnext; k++)
			{
				sum += func(x,fp);
				x += delta;
			}
			t = 0.5*(t + (b-a)*sum/nnext);
			nnext *= 2;

			p = (4*t - tprev)/3;				// the predictor
			if (Math.abs(p-pprev)<TOL*(1+Math.abs(pprev)))
				return p;						// normal exit
			pprev = p;
			tprev = t;
		}
		return -999;						// failed to converge in JMAX iterations
	}
	static double doSimpsonage(double a, double b, double fp[])
	{
		double t, tprev=-999;					
		double p, pprev=-999;					
		int nnext = 1;
		t = 0.5*(b-a)*(funcage(a,fp) + funcage(b,fp));		
		for (int j=1; j<JMAX; j++)
		{
			double delta = (b-a)/nnext;
			double x = a + 0.5*delta;
			double sum = 0.0;					 
			for (int k=0; k<nnext; k++)
			{
				sum += funcage(x,fp);
				x += delta;
			}
			t = 0.5*(t + (b-a)*sum/nnext);
			nnext *= 2;

			p = (4*t - tprev)/3;				
			if (Math.abs(p-pprev)<TOL*(1+Math.abs(pprev)))
				return p;						
			pprev = p;
			tprev = t;
		}
		return -999;						
	}
}
