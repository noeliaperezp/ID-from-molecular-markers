// SNP_ID.c

#include "libhdr"

#define CC 2001   // Maximum N=1000 individuals

int i, nind;

double IDFhatI, IDFhatII, IDFhatIII, IDFhom, IDFexh, IDFroh100, IDFroh1000, IDFroh2000;

struct acc Phe, AFhatI, AFhatII, AFhatIII, AFhom, AFexh, AFroh100, AFroh1000, AFroh2000;

struct covacc PheFhatI, PheFhatII, PheFhatIII, PheFhom, PheFexh, PheFroh100, PheFroh1000, PheFroh2000;
struct covacc FhatIFhatII, FhatIFhatIII, FhatIFhom, FhatIFexh, FhatIFroh100, FhatIFroh1000, FhatIFroh2000;
struct covacc FhatIIFhatIII, FhatIIFhom, FhatIIFexh, FhatIIFroh100, FhatIIFroh1000, FhatIIFroh2000;
struct covacc FhatIIIFhom, FhatIIIFexh, FhatIIIFroh100, FhatIIIFroh1000, FhatIIIFroh2000;
struct covacc FhomFexh, FhomFroh100, FhomFroh1000, FhomFroh2000;
struct covacc FexhFroh100, FexhFroh1000, FexhFroh2000;
struct covacc Froh100Froh1000, Froh100Froh2000;
struct covacc Froh1000Froh2000;

FILE *fFfile, *foutfile;

main()
{
	foutfile = fopen ("outfileID","w");

	getintandskip("NIND :",&nind,10,10000);
	int status = system("bash shell_F_values");
	regression_pm_F();
	printing();
	return(0);
}

/* **************************************************************************** */

regression_pm_F ()
{
	double w, genvalue[CC], FhatI[CC], FhatII[CC], FhatIII[CC], Fhom[CC], Fexh[CC], Froh100[CC], Froh1000[CC], Froh2000[CC];

	fFfile = fopen ("data.F","r");

	for (i=1; i<=nind; i++)
	{
		fscanf(fFfile,"%lf", &w);
		genvalue[i] = w;

		fscanf(fFfile,"%lf", &w);
		FhatI[i] = w;

		fscanf(fFfile,"%lf", &w);
		FhatII[i] = w;

		fscanf(fFfile,"%lf", &w);
		FhatIII[i] = w;

		fscanf(fFfile,"%lf", &w);
		Fhom[i] = w;

		fscanf(fFfile,"%lf", &w);
		Fexh[i] = w;

		fscanf(fFfile,"%lf", &w);
		Froh100[i] = w;

		fscanf(fFfile,"%lf", &w);
		Froh1000[i] = w;

		fscanf(fFfile,"%lf", &w);
		Froh2000[i] = w;
	}

	close(fFfile);

	for (i=1; i<=nind; i++)
	{
		accum (&Phe, genvalue[i]);
		accum (&AFhatI, FhatI[i]);
		accum (&AFhatII, FhatII[i]);
		accum (&AFhatIII, FhatIII[i]);
		accum (&AFhom, Fhom[i]);
		accum (&AFexh, Fexh[i]);
		accum (&AFroh100, Froh100[i]);
		accum (&AFroh1000, Froh1000[i]);
		accum (&AFroh2000, Froh2000[i]);

		covaccum (&PheFhatI, genvalue[i], FhatI[i]);
		covaccum (&PheFhatII, genvalue[i], FhatII[i]);
		covaccum (&PheFhatIII, genvalue[i], FhatIII[i]);
		covaccum (&PheFhom, genvalue[i], Fhom[i]);
		covaccum (&PheFexh, genvalue[i], Fexh[i]);
		covaccum (&PheFroh100, genvalue[i], Froh100[i]);
		covaccum (&PheFroh1000, genvalue[i], Froh1000[i]);
		covaccum (&PheFroh2000, genvalue[i], Froh2000[i]);

		covaccum (&FhatIFhatII, FhatI[i], FhatII[i]);
		covaccum (&FhatIFhatIII, FhatI[i], FhatIII[i]);
		covaccum (&FhatIFhom, FhatI[i], Fhom[i]);
		covaccum (&FhatIFexh, FhatI[i], Fexh[i]);
		covaccum (&FhatIFroh100, FhatI[i], Froh100[i]);
		covaccum (&FhatIFroh1000, FhatI[i], Froh1000[i]);
		covaccum (&FhatIFroh2000, FhatI[i], Froh2000[i]);

		covaccum (&FhatIIFhatIII, FhatII[i], FhatIII[i]);
		covaccum (&FhatIIFhom, FhatII[i], Fhom[i]);
		covaccum (&FhatIIFexh, FhatII[i], Fexh[i]);
		covaccum (&FhatIIFroh100, FhatII[i], Froh100[i]);
		covaccum (&FhatIIFroh1000, FhatII[i], Froh1000[i]);
		covaccum (&FhatIIFroh2000, FhatII[i], Froh2000[i]);

		covaccum (&FhatIIIFhom, FhatIII[i], Fhom[i]);
		covaccum (&FhatIIIFexh, FhatIII[i], Fexh[i]);
		covaccum (&FhatIIIFroh100, FhatIII[i], Froh100[i]);
		covaccum (&FhatIIIFroh1000, FhatIII[i], Froh1000[i]);
		covaccum (&FhatIIIFroh2000, FhatIII[i], Froh2000[i]);

		covaccum (&FhomFexh, Fhom[i], Fexh[i]);
		covaccum (&FhomFroh100, Fhom[i], Froh100[i]);
		covaccum (&FhomFroh1000, Fhom[i], Froh1000[i]);
		covaccum (&FhomFroh2000, Fhom[i], Froh2000[i]);

		covaccum (&FexhFroh100, Fexh[i], Froh100[i]);
		covaccum (&FexhFroh1000, Fexh[i], Froh1000[i]);
		covaccum (&FexhFroh2000, Fexh[i], Froh2000[i]);

		covaccum (&Froh100Froh1000, Froh100[i], Froh1000[i]);
		covaccum (&Froh100Froh2000, Froh100[i], Froh2000[i]);

		covaccum (&Froh1000Froh2000, Froh1000[i], Froh2000[i]);
	}

	IDFhatI = covariance(&PheFhatI) / variance(&AFhatI);
	IDFhatII = covariance(&PheFhatII) / variance(&AFhatII);
	IDFhatIII = covariance(&PheFhatIII) / variance(&AFhatIII);
	IDFhom = covariance(&PheFhom) / variance(&AFhom);
	IDFexh = covariance(&PheFexh) / variance(&AFexh);
	IDFroh100 = covariance(&PheFroh100) / variance(&AFroh100);
	IDFroh1000 = covariance(&PheFroh1000) / variance(&AFroh1000);
	IDFroh2000 = covariance(&PheFroh2000) / variance(&AFroh2000);
}


/* ***************************************************************** */


printing()
{
	fprintf(foutfile, "AFhatI=%6.4f    AFhatII=%6.4f    AFhatIII=%6.4f\n", accmean(&AFhatI), accmean(&AFhatII), accmean(&AFhatIII)); 
	fprintf(foutfile, "AFhom=%6.4f    AFexh=%6.4f\n", accmean(&AFhom), accmean(&AFexh)); 
	fprintf(foutfile, "AFroh100=%6.4f    AFroh1000=%6.4f    AFroh2000=%6.4f\n", accmean(&AFroh100), accmean(&AFroh1000), accmean(&AFroh2000)); 

	fprintf(foutfile, "VFhatI=%6.4f    VFhatII=%6.4f    VFhatIII=%6.4f\n", variance(&AFhatI), variance(&AFhatII), variance(&AFhatIII)); 
	fprintf(foutfile, "VFhom=%6.4f    VFexh=%6.4f\n", variance(&AFhom), variance(&AFexh)); 
	fprintf(foutfile, "VFroh100=%6.4f    VFroh1000=%6.4f    VFroh2000=%6.4f\n", variance(&AFroh100), variance(&AFroh1000), variance(&AFroh2000)); 

	fprintf(foutfile, "            FhatI    FhatII   FhatIII   Fhom    Fexh     Froh100  Froh1000 Froh2000\n");
	fprintf(foutfile, "Phe        %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f\n",
		correlation(&PheFhatI), correlation(&PheFhatII), correlation(&PheFhatIII), correlation(&PheFhom), correlation(&PheFexh), correlation(&PheFroh100), correlation(&PheFroh1000), correlation(&PheFroh2000)); 
	fprintf(foutfile, "FhatI               %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f\n",
		correlation(&FhatIFhatII), correlation(&FhatIFhatIII), correlation(&FhatIFhom), correlation(&FhatIFexh), correlation(&FhatIFroh100), correlation(&FhatIFroh1000), correlation(&FhatIFroh2000)); 
	fprintf(foutfile, "FhatII                        %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f\n",
		correlation(&FhatIIFhatIII), correlation(&FhatIIFhom), correlation(&FhatIIFexh), correlation(&FhatIIFroh100), correlation(&FhatIIFroh1000), correlation(&FhatIIFroh2000)); 
	fprintf(foutfile, "FhatIII                                %6.4f   %6.4f   %6.4f   %6.4f   %6.4f\n",
		correlation(&FhatIIIFhom), correlation(&FhatIIIFexh), correlation(&FhatIIIFroh100), correlation(&FhatIIIFroh1000), correlation(&FhatIIIFroh2000)); 
	fprintf(foutfile, "Fhom                                            %6.4f   %6.4f   %6.4f   %6.4f\n",
		correlation(&FhomFexh), correlation(&FhomFroh100), correlation(&FhomFroh1000), correlation(&FhomFroh2000)); 
	fprintf(foutfile, "Fexh                                                     %6.4f   %6.4f   %6.4f\n",
		correlation(&FexhFroh100), correlation(&FexhFroh1000), correlation(&FexhFroh2000)); 
	fprintf(foutfile, "Froh100                                                           %6.4f   %6.4f\n", correlation(&Froh100Froh1000), correlation(&Froh100Froh2000)); 
	fprintf(foutfile, "Froh1000                                                                   %6.4f\n", correlation(&Froh1000Froh2000));

	fprintf(foutfile, "IDFhatI=%6.4f    IDFhatII=%6.4f    IDFhatIII=%6.4f\n", IDFhatI, IDFhatII, IDFhatIII); 
	fprintf(foutfile, "IDFhom=%6.4f    IDFexh=%6.4f\n", IDFhom, IDFexh); 
	fprintf(foutfile, "IDFroh100=%6.4f    IDFroh1000=%6.4f    IDFroh2000=%6.4f\n", IDFroh100, IDFroh1000, IDFroh2000); 

	return(0);
}

/* ********************************************************************************************* */
