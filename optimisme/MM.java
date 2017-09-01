package optimisme;

import java.util.ArrayList;
import java.util.Arrays;


import optimisme.PSFPreparator;
import optimisme.MMCal;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.measure.Calibration;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import edu.emory.mathcs.jtransforms.fft.DoubleFFT_3D;
import org.ejml.simple.SimpleMatrix;

/**
 * Title:        OPTIMIST --- A 3D  implementation for ImageJ
 * Description: Optimism is a deconvolution method based on the Majorization
 * minimization algorithm, devoted to bi-photon' microscopic images.
 * 
 * Copyright:    Copyright (c) 2017
 * Company:      Labex Archimède, I2M
 *               AMU Aix Marseille Université
 * @author       Dominique Benielli
 * @author       Caroline Chaux-Moulin
 * @author       Sandrine Anthoine
 * @version 1.2
 */
public class MM 
{
	// 
	private double[][][] dataIm = null;
	private double[][][] dataImDeconvol = null;
	private double[][][] dataPsfInit;
	private double[][][] dataPsf = null;
	private MMParameter mmparam;
	private int[] dataDimIm = new int[3];
	private boolean function_handle = false;
	

	
	/**
	 * This Class instantiate the Majorization Minimization
	 * Algorithms, populate fields dataPsf with the resampleing PSF
	 * and dataImDeconvol and the decovolued image, result of the alrorithms
	 */
	public MM(){

	}

	/**
	 * Constructor of MM Class
	 * @param sourceIm ImageStack input image
	 * @param sourcePSF ImageStack input PSF image
	 * @param mmparam MMParameter contains parameters values for the process
	 */
	public MM( ImageStack sourceIm, ImageStack sourcePSF, MMParameter mmparam ){
		IJ.showProgress(0.0);
		this.process(sourceIm, sourcePSF, mmparam);
		IJ.showProgress(1.0);
	}

	/**
	 * Get the dataIm 3D data Image
	 * @return 3D data image
	 */
	public double[][][] getDataIm() {
		return dataIm;
	}
	
	/**
	 * Get the dataImDeconvol 3D data Image processed 
	 * @return 3D data image, processed
	 */
	public double[][][] getDataImDeconvol() {
		return dataImDeconvol;
	}

	/**
	 * Set the dataIm 3D data Image
	 * @param dataIm 3D data image
	 */
	public void setDataIm(double[][][] dataIm) {
		this.dataIm = dataIm;
	}

	/**
	 * Get the dataPsfInit 3D data PSF input Image
	 * @return 3D data PSF input image
	 */
	public double[][][] getDataPsfInit() {
		return dataPsfInit;
	}

	/**
	 * Set the dataPsfInit 3D data PSF input Image
	 * @param dataPsfInit data PSF input image
	 */
	public void setDataPsfInit(double[][][] dataPsfInit) {
		this.dataPsfInit = dataPsfInit;
	}

	/**
	 * Get the dataPsf  3D data PSF Image
	 * @return dataPsf data PSF image processed
	 */
	public double[][][] getDataPsf() {
		return dataPsf;
	}

	/**
	 * Set the dataPsf  3D data PSF Image
	 * @param dataPsf
	 */
	public void setDataPsf(double[][][] dataPsf) {
		this.dataPsf = dataPsf;
	}

	/**
	 * Get the MMParameter values
	 * @return MMParameter parameters of algorithms
	 */
	public MMParameter getMmparam() {
		return mmparam;
	}

	/**
	 * Set the MMParameter values
	 * @param mmparam MMParameter of algorithm
	 */
	public void setMmparam(MMParameter mmparam) {
		this.mmparam = mmparam;
	}

	/**
	 * Get image data Table  containing the 3 dimensions of image
	 * (1D of three integers)
	 * @return dataDimIm table of image data
	 */
	public int[] getDataDimIm() {
		return dataDimIm;
	}

	/**
	 * Set image data Table  containing the 3 dimensions of image
	 * @param dataDimIm table of 1D of three integers number
	 */
	public void setDataDimIm(int[] dataDimIm) {
		this.dataDimIm = dataDimIm;
	}


	/**
	 * This method process the MM algorithm
	 * @param sourceIm ImageStack image input
	 * @param sourcePSF ImageStack PSF input image
	 * @param mmparam  MMParameter parameters of the algorithm
	 */
	public void process( ImageStack sourceIm, ImageStack sourcePSF, MMParameter mmparam )
	{
		this.dataDimIm[0] =  (int) sourceIm.getWidth();
		this.dataDimIm[1] = sourceIm.getHeight();
		this.dataDimIm[2] = sourceIm.getSize();
		this.mmparam = mmparam;
		double[][][] newdataPSF = null;    
		this.dataPsfInit =  populateData(sourcePSF); 
		this.mmparam = mmparam;

		int[] sizeOUT = {Integer.MAX_VALUE, Integer.MAX_VALUE, Integer.MAX_VALUE};

		/* ****************************
		 * Boundary Initialization 
		 * ****************************/
		 int xmin = 0;
		 int xmax = Integer.MAX_VALUE;

		 /* ******************************************
		  * PSF preparator if necessary
		  *******************************************/ 
		 IJ.showProgress(0.05);
		 if (!this.compareResolution(mmparam.resIM,mmparam.resPSF)) {
			 PSFPreparator psfPreparator = new PSFPreparator();      	
			 this.dataPsf = psfPreparator.process(dataPsfInit, mmparam, sizeOUT );
			 this.plotCalculatedPsf(); 
			 IJ.showProgress(0.49);
		 }
		 else {
			 double normpsf = MMCal.sum( MMCal.concatAll(this.dataPsfInit))[0];
			 this.dataPsf = MMCal.rDivide (this.dataPsfInit, normpsf );     
		 }
		 /* ***************************************
		  * Image populate
		  *****************************************/
		 this.dataPsfInit = null;
		 this.dataIm =  populateData(sourceIm);   

		 //****************************************
		 // * Normanization of Image data
		 //***************************************       
		 double normH = MMCal.sum( MMCal.concatAll(this.dataIm))[0];
		 this.dataIm = MMCal.rDivide (this.dataIm, normH ); 

		 //***********************************************
		 //* Make a copy for algorithm generalarization 
		 //************************************************
		 double[][][]  X0 = MMCal.copy(this.dataIm);
		 double[][][]  y = this.dataIm;
		 IJ.showProgress(0.5);
                 System.gc();
		 ArrayList l = this.MMAlgorithm(y, X0, this.dataIm, xmin, xmax);
                 System.gc();
         this.dataImDeconvol = (double[][][]) l.get(2);
         this.plotCalculatedImage();
		 IJ.showProgress(1.0);

	}

	//*************************************************
	//* new Psf visualisation and calibration populate
	//*************************************************
	/**
	 * This method realizes the plot of the dataPsf data 
	 * and populates the new image resolution with the same as image
	 */
	private void plotCalculatedPsf(){
		ImageStack psfStack = new ImageStack(this.dataPsf[0][0].length,
				this.dataPsf[0].length);

		for (int z=0 ; z < this.dataPsf.length ; z++)
		{  
			float[][] data = MMCal.converttoFloat2dandReverseDim(dataPsf[z]);
			ImageProcessor ip = new FloatProcessor(data);
			psfStack.addSlice("",ip);
		}

		ImagePlus image = new ImagePlus("Re-sampled PSF", psfStack);			
		Calibration calib = image.getCalibration();
		calib.pixelDepth = this.mmparam.resIM[2];
		calib.pixelHeight = this.mmparam.resIM[1];
		calib.pixelWidth = this.mmparam.resIM[0];
		calib.setUnit(this.mmparam.unitRes);
		image.setCalibration(calib);
		image.show();
	}

	//********************************************
	//* new Image deconvolved visualisation 
	//********************************************
	/**
	 * This method realizes the plot of the dataImDeco data 
	 * and populates the new image resolution with the same as image
	 */
	private void plotCalculatedImage(){
		ImageStack imageStack = new ImageStack(this.dataImDeconvol[0][0].length,
				this.dataImDeconvol[0].length);

		for (int z=0 ; z < this.dataImDeconvol.length ; z++)
		{  
			float[][] data = MMCal.converttoFloat2dandReverseDim(dataImDeconvol[z]);
			ImageProcessor ip = new FloatProcessor(data);
			//ImageProcessor ip = new FloatProcessor(dataPsf[0][0].length ,dataPsf[0].length);
			//	for(int y=0 ; y< dataPsf[z].length ; y++)
			//	{   
			//		float[] data = MMCal.converttoFloat1d(dataPsf[z]);
			//		ip.putRow(0, y, data, dataPsf[z][y].length);
			//	}
			imageStack.addSlice("",ip);
		}

		ImagePlus image = new ImagePlus("Deconvolved image", imageStack);			
		Calibration calib = image.getCalibration();
		calib.pixelDepth = this.mmparam.resIM[2];
		calib.pixelHeight = this.mmparam.resIM[1];
		calib.pixelWidth = this.mmparam.resIM[0];
		calib.setUnit(this.mmparam.unitRes);
		image.setCalibration(calib);
		image.show();
	}
	
	//********************************************************************
	//* compare resolution of image and psf
	//********************************************************************
	/**
	 * This method compares the resolutions of the both resolution of data
	 * if return True the dimensions are the same False otherwise
	 * @param resIm  array of 3 with values of Image resolution in three directions
	 * @param resPSF array of 3 with values of PSF Image resolution in three directions
	 * @return boolean True or False
	 */
	protected boolean compareResolution(double[] resIm, double[] resPSF){

		if ((resIm[0] == resPSF[0]) && 
				(resIm[1] == resPSF[1]) &&
				(resIm[2] == resPSF[2]) )
			return true;
		else return false;
	}

	//********************************************************************
	//* compare dimensions of image and psf
	//********************************************************************   
	/**
	 * This method compares the dimensions of the both input data
	 * if return True the dimensions are the same False otherwise
	 * @param dataIm 3D input data
	 * @param dataPSF 3D input data
	 * @return boolean True or False
	 */
	protected boolean compareDimension(double[][][]  dataIm, double[][][]  dataPSF) {   	
		if (     ( ( dataIm.length != dataPSF.length ) ||
				( dataIm[ 0 ].length != dataPSF[ 0 ].length ) ||
				( dataIm[ 0 ][ 0 ].length != dataPSF[ 0 ][ 0 ].length ) ) )
			return true;

		else
			return false;
	}

	/**
	 * This method setData sets the Inputs ImageStack to a 3D data attributes dataIm and dataPsf 
	 * @param sourceIm ImageStack input Image
	 * @param sourcePSF ImageStack input PSF Image 
	 */
	public void setData( ImageStack sourceIm ,ImageStack sourcePSF)
	{
		if ( sourceIm == null )
			throw new IllegalArgumentException( "Source stack with image cannot be 'null'." );
		if ( sourcePSF == null )
			throw new IllegalArgumentException( "Source stack with PSF cannot be 'null'." );
		this.dataIm = populateData(sourceIm);
		this.dataPsf = populateData(sourcePSF);
	}



	/**
	 * This method populateData converts the Input ImageStack to a 3D data
	 * @param source ImageStack
	 * @return 3D output data
	 */
	public double[][][] populateData(ImageStack source){
		int dataDimX = source.getWidth();
		int dataDimY = source.getHeight();
		int dataDimZ = source.getSize();
		double[][][] data = new double[dataDimZ][dataDimY][dataDimX];

		for ( int z = 0; z < dataDimZ; z++ )
		{
			ImageProcessor imageIm = source.getProcessor( z + 1 );
			if ( source != null )
			{
				for ( int x = 0; x < dataDimX; x++ )
					for ( int y = 0; y < dataDimY; y++ ){
						data[ z ][ y ][ x ] = imageIm.getPixelValue( x, y );
					}
			}
		}
		return data;
	}


	/**
	 * This method performs the main algorithm of Majorization Minimization	
	 * The method return un list of arrays End SignalNoiseRatio, SignalNoiseRatio
	 * 3D output image (x), internal algorithm Critere and NGradiant.
	 * y , x , xbar pointed to the same data, do copy for generalized algorithm
	 * @param y 3D input data of image
	 * @param x input data of image
	 * @param xbar input data of image
	 * @param xmin minimum of x
	 * @param xmax maximum of x
	 * @return ArrayList contains : SNRend, SNR, x, Crit ,NGrad
	 */
	//****************************************************************************
	//* y , x , xbar pointed to the same data, do copy for generalized algorithm
	//*****************************************************************************
	
	@SuppressWarnings("rawtypes")
	public ArrayList MMAlgorithm( double[][][] y, double[][][] x,double[][][] xbar,
			int xmin, int xmax)
	{ 
		IJ.showProgress(0.51);
		//****************************************************************************
		//*INITIALIZATION OF VARIABLES
		//*****************************************************************************
		int NbIt = this.mmparam.NbIt;
		ArrayList lo = new ArrayList();
		int Nx = x[0][0].length;
		int Ny = x[0].length;
		int Nz = x.length;
		double stop = 1E-2;
		int phiXY = this.mmparam.phixy; 	
		int phiZ = this.mmparam.phiz;
		double[] NGrad = new double[NbIt+1];
		double[] Crit = new double[NbIt+1];
		double[] Time = new double[NbIt+1];
		double[] SNR = new double[NbIt+1];
		double[] Err = new double[NbIt+1];

		//****************************************************************************
		//*INITIALIZATION OF ALGORITHM
		//*****************************************************************************    	
		ArrayList li =  this.critere(x, y, xmin, xmax);
	//	Crit[0] = Double.valueOf(li.get(0).toString()); 
		Crit[0] =((double[]) li.get(0))[0];
	//	System.out.print("crit"+String.valueOf(Crit[0])+ "  ");
		double[][][] Grad ;
		Grad = (double[][][]) li.get(1);
		
	//	System.out.print("Grad "+toString(Grad[0][0])+ "  ");
		double[][][] wXY_Vx;
		wXY_Vx = (double[][][]) li.get(2);
		double[][][] wZ_Vx = (double[][][]) li.get(3);
		// Initial criterion value 
		Time[0] = 0;
		NGrad[0] = MMCal.norm(  MMCal.concatAll(Grad));
	//	System.out.print("NGrad "+String.valueOf(NGrad[0])+ "  ");
		if (xbar == null)
			SNR[0] = Double.POSITIVE_INFINITY;
		else {

			Err[0] = MMCal.sum( MMCal.square(MMCal.minus( MMCal.concatAll(x),MMCal.concatAll(xbar))) )[0]
					/ MMCal.sum(MMCal.square( MMCal.concatAll(xbar) ) )[0];    	      	    
			SNR[0] = 10*Math.log10(1/Err[0]);
		}
		double SNRend = Double.POSITIVE_INFINITY;
		double[][][] dx = null;
		double[][][] Hdx = null;
		double[][][] Vvdx = null;
		double[][][] Vhdx = null;
		double[][][] Vtdx = null;
		double B;
		double timeStep = ((1.0-0.52)/(double)(NbIt));
		System.gc();
		IJ.showProgress(0.52);
		for (int k = 0 ; k< NbIt ; k++){
			
			//IJ.showProgress(0.52+(double)k*timeStep);
			//tic

			// Stopping test
			if(NGrad[k] < stop){
				break;
			}
			if (k == Math.floorDiv(NbIt, 9) || NGrad[k]/stop < 10-5) {
				IJ.showProgress(0.6);
			}
			if (k == Math.floorDiv(NbIt, 8)|| NGrad[k]/stop < 10-4) {
				IJ.showProgress(0.65);
			}
			if (k == Math.floorDiv(NbIt, 7) || NGrad[k]/stop < 10-3) {
				IJ.showProgress(0.7);
			}			
			if (k == Math.floorDiv(NbIt, 6)|| NGrad[k]/stop < 10-2) {
				IJ.showProgress(0.75);
			}
			if (k == Math.floorDiv(NbIt, 5)|| NGrad[k]/stop < 10-1) {
				IJ.showProgress(0.8);
			}
			if (k == Math.floorDiv(NbIt, 4)) {
				IJ.showProgress(0.85);
			}
			if (k == Math.floorDiv(NbIt, 3)) {
				IJ.showProgress(0.9);
			}	
			if (k == Math.floorDiv(NbIt, 2)) {
				IJ.showProgress(0.95);
			}			
			//IJ.showProgress(0.52+(double)k*timeStep);

			ArrayList loo = MMCal.Voperator(Grad);
		//	ArrayList loo = MMCal.Voperator(this.dataPsf);
			double[][][] Vvg = (double[][][]) loo.get(1);
			double[][][] Vhg = (double[][][]) loo.get(0);
			double[][][] Vtg = (double[][][]) loo.get(2);
				
	//		System.out.print("Vvg "+toString(Vvg));

			double[][][] Hg = this.apply_PSFvar3D(Grad);
			if (k==0){

				// %Dir   = -Grad;         
				B = this.majorante1(x, wXY_Vx, wZ_Vx, MMCal.dotMult(MMCal.copy(Grad),-1.0 ), Vvg, Vhg, Vtg, Hg, xmin,xmax);
				double s = MMCal.sum(MMCal.square( MMCal.concatAll(Grad) ))[0]/B;
				dx = MMCal.dotMult(MMCal.copy(Grad), -s);
				Hdx = MMCal.dotMult(MMCal.copy(Hg), -s);
				Vvdx = MMCal.dotMult(MMCal.copy(Vvg),-s);

				Vhdx = MMCal.dotMult(MMCal.copy(Vhg),-s);
				Vtdx = MMCal.dotMult(MMCal.copy(Vtg),-s);
			}

			else
			{

				double[][][] Dir1 = MMCal.dotMult(MMCal.copy(Grad),-1.0 );
				double[][][] Dir2 = dx;

				double[][][] Hd1 = MMCal.dotMult(MMCal.copy(Hg), -1.0); 
				double[][][] Hd2 = MMCal.copy(Hdx);
				double[][][] Vvd1 =  MMCal.dotMult(MMCal.copy(Vvg),-1.0);
				double[][][] Vvd2 = MMCal.copy(Vvdx);
				double[][][] Vhd1 =  MMCal.dotMult(MMCal.copy(Vhg),-1.0); 
				double[][][] Vhd2 = MMCal.copy(Vhdx);
				double[][][] Vtd1 =  MMCal.dotMult(MMCal.copy(Vtg),-1.0);
				double[][][] Vtd2 = MMCal.copy(Vtdx);

				double[][] B2 = this.majorante2(x,wXY_Vx,wZ_Vx,Dir1,Dir2,Vvd1,Vhd1,Vtd1,Vvd2,Vhd2,Vtd2,Hd1,Hd2,xmin,xmax);
			//	System.out.print("B2 "+ this.toString(B2)+ "  ");
				double[] concatGrad = MMCal.concatAll(Grad);
				double d1 = MMCal.sum(   MMCal.dotMult(MMCal.concatAll(Dir1),concatGrad) )[0];
				double d2 = MMCal.sum(   MMCal.dotMult(MMCal.concatAll(Dir2),concatGrad) )[0];

				SimpleMatrix a = new SimpleMatrix(B2);

				double[][] d =  new double[2][1];
				d[0][0] = d1;
				d[1][0] = d2;               
				SimpleMatrix dmat = new SimpleMatrix(d);
				SimpleMatrix i =  a.pseudoInverse().mult(dmat).scale(-1);

				dx = MMCal.add( MMCal.dotMult(Dir1, i.get(0, 0)) , MMCal.dotMult(Dir2, i.get(1, 0))); // Dir1 Dir2 compromised

				Hdx =  MMCal.add( MMCal.dotMult(Hd1, i.get(0, 0)) , MMCal.dotMult(Hd2, i.get(1, 0))); // Hd1 Hd2 compromised
				Vvdx =  MMCal.add( MMCal.dotMult(Vvd1 , i.get(0, 0)) , MMCal.dotMult(Vvd2, i.get(1, 0))); // Vvd1 Vvd2 compromised
				Vhdx = MMCal.add( MMCal.dotMult(Vhd1, i.get(0, 0)) , MMCal.dotMult(Vhd2, i.get(1, 0))) ; // Vhd1 Vhd2 compromised
				Vtdx =  MMCal.add( MMCal.dotMult(Vtd1, i.get(0, 0)) , MMCal.dotMult(Vtd2, i.get(1, 0)));   	//  Vtd1   Vtd2 compromised   
			}

			// update
			x = MMCal.add(x, dx);

			li = critere(x, y, xmin, xmax);
			Crit[k+1] = ((double[])  li.get(0))[0]; 	
	//		System.out.print("crit "+String.valueOf(Crit[k+1])+ "  ");
			Grad = (double[][][]) li.get(1);
			wXY_Vx = (double[][][]) li.get(2);
			wZ_Vx = (double[][][]) li.get(3);

			//   Time[k+1] = toc;
			if(xbar!=null){
				Err[k+1] = MMCal.sum( MMCal.square(MMCal.minus( MMCal.concatAll(MMCal.copy(x)),MMCal.concatAll(xbar))) )[0]
						/ MMCal.sum(MMCal.square( MMCal.concatAll(xbar) ) )[0];  
				SNR[k+1] =  10*Math.log10(1/Err[k+1]); 

			}
			NGrad[k+1] = MMCal.norm(MMCal.concatAll(Grad));

		}
		SNRend = SNR[SNR.length-1];

		lo.add(SNRend);
		lo.add(SNR);
		lo.add(x);
		lo.add(Crit);	
		lo.add(NGrad);
		return lo;
	}

	
	/**
	 * This method apply_PSFvar3D performs internal computation 
	 * Implementation only for on PSF
	 * @param x 3D data image
	 * @return 3D data  of Fourier convolution
	 */
	public double[][][] apply_PSFvar3D(double[][][] x){
		double[][][] result = null;
		if (function_handle){
			;
		}
		else {
			result = MMCal.conv3D_fourier(x, this.dataPsf);
		}

		return result;
	}

	/**
	 * This method  apply_PSFadjvar3D performs internal computation 
	 * Implementation only for on PSF
	 * @param x 3D data image
	 * @return 3D data  of PSFadjvar3D
	 */
	public double[][][] apply_PSFadjvar3D(double[][][] x){
		double[][][] result = null;
		double[][][] a = null;
		int N1 = x[0][0].length;
		int N2 = x[0].length;
		int N3 = x.length;
		if (function_handle){
			;
		}
		else {
			a = this.dataPsf;

		}

		if (!function_handle){

			int[] p = new int[3];
			p[0] = Math.floorDiv(a[0][0].length,2);
			p[1] = Math.floorDiv(a[0].length,2);
			p[2] = Math.floorDiv(a.length,2);
			double[][][] Dext = new double[N3+2*p[2]] [N2+2*p[1]] [N1+2*p[0]];

			for (int k = 0; k < N3  ; k++)
				for (int j = 0; j < N2 ; j++)
					for (int i = 0 ; i < N1  ; i++)
						Dext[k+p[2]][j+p[1]][i+p[0]] = x[k][j][i];
			double[][][] Hcomplex = MMCal.conj(MMCal.freqfilt3D(a, 2*p[0]+N1, 2*p[1]+N2,  2*p[2]+N3));


			DoubleFFT_3D FFTCal = new DoubleFFT_3D(Dext.length,Dext[0].length,Dext[0][0].length);
			double[][][] DextComplex = MMCal.complexCopy(Dext);
			FFTCal.complexForward(DextComplex);

			DextComplex = MMCal.dotMultComplex(DextComplex, Hcomplex);
			FFTCal.complexInverse(DextComplex, true);


			int[] indexX = new int[N1];
			int[] indexY = new int[N2];
			int[] indexZ = new int[N3];
			for(int i = 0 ; i < N1 ; i++)
				indexX[i] = p[0] +  i;
			for(int i = 0 ; i < N2 ; i++)
				indexY[i] = p[1] +  i;
			for(int i = 0 ; i < N3 ; i++)
				indexZ[i] = p[2] +  i;
			double[][][] HstarxDouble = MMCal.toDoubleReal(DextComplex);

			DextComplex = null;
			HstarxDouble = MMCal.getIndex(HstarxDouble, indexX, indexY, indexZ);
			return HstarxDouble;

		}
		else{
			double p3 = (a.length-1.0)/2.0;
			//for ( int z = 0 ; z < x.length ; z++ )
				//	{   
				//	for (int n3 = Math.max(0, (int) (z-p3))  ; n3 < Math.min(x.length, (int)(sa3+z-1-p3)) ; n3++){
					//		double a = MMCal.concatAll(this.dataPsf)[n3];
			//		bF2 = conv2Dadjoint_fourier(x(:,:,n3),a(:,:,n3-z+p3+1));
			//		}
			//	}
			return result;
		}
		
	}

	/// attention a finir instantiation et input a verifier
	/**
	 * This method   majorante1 evaluates the majorantes criteria 1
	 * @param x 3D data image
	 * @param wXY_Vx 3D data gradient
	 * @param wZ_Vx 3D data gradient
	 * @param d1 3D data 
	 * @param Vvd1 3D data 
	 * @param Vhd1 3D data 
	 * @param Vtd1 3D data 
	 * @param Hd1 3D data 
	 * @param xmin Integer parameter for min
	 * @param xmax Integer parameter for max
	 * @return double evaluation of majoration
	 */
	public double majorante1(double[][][] x,double[][][] wXY_Vx,double[][][] wZ_Vx,
			double[][][] d1,double[][][] Vvd1,double[][][] Vhd1,double[][][] Vtd1,
			double[][][] Hd1,int xmin,int xmax){
		double lambdaXY = this.mmparam.regul;
		double lambdaZ = this.mmparam.regulZ;
		double DthtHD = MMCal.sum( MMCal.square( MMCal.concatAll(Hd1) ) )[0];
	//	double[][][] lDtVtWVD =  new double[Vvd1.length][Vvd1[0].length][Vvd1[0][0].length];
		double lDtVtWVD ;
		double DtVtWVD = 0;
		for (int k = 0; k < Vvd1.length  ; k++)
			for (int j = 0 ; j < Vvd1[0].length ; j++)
				for ( int i = 0 ; i < Vvd1[0][0].length  ; i++){
					lDtVtWVD = lambdaXY * (Vvd1[k][j][i] * (wXY_Vx[k][j][i] * Vvd1[k][j][i]) +
							Vhd1[k][j][i] *( wXY_Vx[k][j][i] *Vhd1[k][j][i]))
							+ lambdaZ * (Vtd1[k][j][i] * (wZ_Vx[k][j][i] *Vtd1[k][j][i]));
					DtVtWVD += lDtVtWVD;
				}


		//double DtVtWVD = MMCal.sum(MMCal.concatAll(lDtVtWVD))[0];

		double[] d1min = MMCal.getWithFilter(d1, MMCal.lte(xmin, x) );  

		double[] d1max =  MMCal.getWithFilter(d1, MMCal.gte(xmax, x) );

		double B = DthtHD + DtVtWVD + 
				this.mmparam.eta * MMCal.sum( MMCal.square( d1min ) )[0]
						+ MMCal.sum(  d1max)[0];   // d1min corrupted
		//System.out.print("B "+String.valueOf(B )+ "  ");
		return B;
	}

	/// attention a finir instantiation et input a verifier

	/**	 
	 * This method   majorante2 evaluates the majorantes criteria 2
	 * @param x 3D data image
	 * @param wXY_Vx 3D data gradient
	 * @param wZ_Vx 3D data gradient
	 * @param d1 3D data 
	 * @param d2 3D data 
	 * @param Vvd1 3D data 
	 * @param Vhd1 3D data 
	 * @param Vtd1 3D data 
	 * @param Vvd2 3D data 
	 * @param Vhd2 3D data 
	 * @param Vtd2 3D data 
	 * @param Hd1 3D data 
	 * @param Hd2 3D data 
	 * @param xmin Integer parameter for min
	 * @param xmax Integer parameter for max
	 * @return 2by2 Matrix of number evaluation of majoration 2
	 */
	public double[][] majorante2(double[][][] x,double[][][] wXY_Vx,double[][][] wZ_Vx,
			double[][][] d1,double[][][] d2,double[][][] Vvd1,double[][][] Vhd1,double[][][] Vtd1,
			double[][][] Vvd2,double[][][] Vhd2,double[][][] Vtd2,double[][][] Hd1,double[][][] Hd2,
			double xmin, double xmax){

		double lambdaXY = this.mmparam.regul;
		double lambdaZ = this.mmparam.regulZ;

		//% d1td1 = sum(d1(:).^2);
		//% d2td2 = sum(d2(:).^2);
		//% d1td2 = sum(d1(:).*d2(:));

		double[] interHd1 = MMCal.concatAll(Hd1);
		double[] interHd2 = MMCal.concatAll(Hd2);
		double d1thtHd1 = MMCal.sum( MMCal.square( Arrays.copyOf(interHd1,interHd1.length)) )[0];
		double d2thtHd2 = MMCal.sum( MMCal.square(Arrays.copyOf(interHd2,interHd2.length)))[0];   	
		double d1thtHd2 = MMCal.sum(MMCal.dotMult(interHd1, interHd2))[0]; // interHd1 compromised
	//	 System.out.print("d1thtHd1 "+String.valueOf(d1thtHd1 )+ "  ");
	//     System.out.print("d2thtHd2 "+String.valueOf(d2thtHd2 )+ "  ");
	//     System.out.print("d1thtHd2 "+String.valueOf(d1thtHd2)+ " \n ");
	//	double[][][] ld1tVtWVd1 = new double[ Vvd1.length][ Vvd1[0].length][ Vvd1[0][0].length];
	//	double[][][] ld1tVtWVd2 = new double[ Vvd1.length][ Vvd1[0].length][ Vvd1[0][0].length];
	//	double[][][] ld2tVtWVd2 = new double[ Vvd1.length][ Vvd1[0].length][ Vvd1[0][0].length];
	    double ld1tVtWVd1;
	    double ld1tVtWVd2;
	    double ld2tVtWVd2;
		double d1tVtWVd1 = 0;
		double d1tVtWVd2 = 0;
		double d2tVtWVd2 = 0;

		
		for (int k = 0; k < Vvd1.length  ; k++)
			for (int j = 0 ; j < Vvd1[0].length ; j++)
				for ( int i = 0 ; i < Vvd1[0][0].length  ; i++){
					ld1tVtWVd1 = lambdaXY *( Vvd1[k][j][i] * (wXY_Vx[k][j][i] * Vvd1[k][j][i])
							+ Vhd1[k][j][i] *(wXY_Vx[k][j][i] * Vhd1[k][j][i]))
							+ lambdaZ *  (Vtd1[k][j][i] * (wZ_Vx[k][j][i] * Vtd1[k][j][i])); 				
					d1tVtWVd1 += ld1tVtWVd1;
					ld1tVtWVd2= lambdaXY * (Vvd1[k][j][i]*(wXY_Vx[k][j][i]*Vvd2[k][j][i])
							+Vhd1[k][j][i] *(wXY_Vx[k][j][i] * Vhd2[k][j][i]))
							+ lambdaZ * (Vtd1[k][j][i] * (wZ_Vx[k][j][i] * Vtd2[k][j][i]));
					d1tVtWVd2 += ld1tVtWVd2;
					ld2tVtWVd2 = lambdaXY * (Vvd2[k][j][i] * (wXY_Vx[k][j][i] * Vvd2[k][j][i] ) 
							+ Vhd2[k][j][i] * (wXY_Vx[k][j][i] * Vhd2[k][j][i]))
							+ lambdaZ * (Vtd2[k][j][i] * (wZ_Vx[k][j][i] * Vtd2[k][j][i]));
					 d2tVtWVd2 += ld2tVtWVd2;
				}

	//	double d1tVtWVd1 = MMCal.sum(MMCal.concatAll(ld1tVtWVd1))[0];
	//	double d1tVtWVd2 = MMCal.sum(MMCal.concatAll(ld1tVtWVd2))[0];
	//	double d2tVtWVd2 = MMCal.sum(MMCal.concatAll(ld2tVtWVd2))[0];
 // /      System.out.print( "d1tVtWVd1 "+String.valueOf(d1tVtWVd1)+ "  ");
 ///       System.out.print(" d1tVtWVd2 " +String.valueOf(d1tVtWVd2)+ "  ");
  ///      System.out.print("d2tVtWVd2 "+String.valueOf(d2tVtWVd2)+ " \n ");
		//%d1UNmind1 = d1.*(UNmin.*d1); d1UNmaxd1 = d1.*(UNmax.*d1); 

		double[] d1min = MMCal.getWithFilter(d1, MMCal.lte(xmin, x) ); 
		double[] d2min = MMCal.getWithFilter(d2, MMCal.lte(xmin, x) ); 
		double[] d1max =  MMCal.getWithFilter(d1, MMCal.gte(xmax, x) );
		double[] d2max =  MMCal.getWithFilter(d2, MMCal.gte(xmax, x) );

		double[] d1UNmind1 = MMCal.square( Arrays.copyOf(d1min,d1min.length ) );  
		double[] d1UNmaxd1 = MMCal.square( Arrays.copyOf(d1max,d1max.length));

		double[] d2UNmind2 = MMCal.square( Arrays.copyOf(d2min,d2min.length)) ; 
		double[] d2UNmaxd2 = MMCal.square( Arrays.copyOf(d2max,d2max.length) ); 

		double[] d1UNmind2 = MMCal.dotMult( d1min, d2min); // d1min compromised
		double[] d1UNmaxd2 = MMCal.dotMult( d1max, d2max); // d1max compromised
		//% B11 = d1thtHd1 + d1tVtWVd1  + 2*eta.*d1td1;
		//% B12 = d1thtHd2 + d1tVtWVd2  + 2*eta.*d1td2;
		//% B22 = d2thtHd2 + d2tVtWVd2  + 2*eta.*d2td2;

		double B11 = d1thtHd1 + d1tVtWVd1  + 
				this.mmparam.eta * MMCal.sum(d1UNmind1)[0] + 
				MMCal.sum(d1UNmaxd1)[0];
		double B12 = d1thtHd2 + d1tVtWVd2  + 
				this.mmparam.eta * MMCal.sum(d1UNmind2)[0] + 
				MMCal.sum(d1UNmaxd2)[0];
		double B22 = d2thtHd2 + d2tVtWVd2  + 
				this.mmparam.eta * MMCal.sum(d2UNmind2)[0] + 
				MMCal.sum(d2UNmaxd2)[0];
		double[][] B = new double[2][2];
		B[0][0] = B11;
		B[0][1] = B[1][0] = B12;
		B[1][1] = B22;
	//	 System.out.print(String.valueOf(B11)+ "  ");
	//	 System.out.print(String.valueOf(B12)+ "  ");
	//	 System.out.print(String.valueOf(B22)+ "  ");
		return B;
	}   

	/**	 
	 * This method  critere  evaluates the criteria of algorithm
	 * @param x 3D data image
	 * @param y 3D data image
     * @param xmin Integer parameter for min
	 * @param xmax Integer parameter for max
	 * @return 	Array list containing;  F (number) , dF : 3D data, wXY_Vx : 3D data , wZ_Vx 3D : data
	 */
	public ArrayList critere(double[][][] x, double[][][] y, int xmin, int xmax){
		ArrayList li = new ArrayList();
		double[][][] Hx = apply_PSFvar3D(x);
		int phiXY = this.mmparam.phixy;

		int phiZ = this.mmparam.phiz;
		double lambdaXY = this.mmparam.regul;
		double lambdaZ = this.mmparam.regulZ;
		double deltaXY = this.mmparam.T;
		double deltaZ = this.mmparam.TZ;
		ArrayList liV = MMCal.Voperator(x);
		double[][][] Vvx = (double[][][]) liV.get(1);
		double[][][] Vhx = (double[][][]) liV.get(0);
		double[][][] Vtx = (double[][][]) liV.get(2);		
		
		double[][][] Vx =  MMCal.sqrt( MMCal.add(
				MMCal.square(MMCal.copy(Vvx)), MMCal.square(MMCal.copy(Vhx)) ));
		double[][][] phiXY_Vx = null;
		double[][][] wXY_Vx = null;
		double[][][] phiZ_Vx = null;
		double[][][] wZ_Vx = null;
		double[][][] sqVx =  MMCal.square( MMCal.copy(Vx));
		double[][][] sqVtx =  MMCal.square( MMCal.copy(Vtx));
		switch (phiXY) {
		case 1:
			phiXY_Vx = MMCal.phi_sqV_1(sqVx, deltaXY);      	 
			wXY_Vx = MMCal.w_V_1(sqVx, deltaXY);       			       
			break;
		case 2:
			phiXY_Vx = MMCal.phi_sqV_2(sqVx, deltaXY);      	 
			wXY_Vx = MMCal.w_V_2(sqVx, deltaXY);
			break;
		case 3:
			phiXY_Vx = MMCal.phi_sqV_3(sqVx, deltaXY);      	 
			wXY_Vx = MMCal.w_V_3(sqVx, deltaXY);
			break;
		case 4:
			phiXY_Vx = MMCal.phi_sqV_4(sqVx, deltaXY);      	 
			wXY_Vx = MMCal.w_V_4(sqVx, deltaXY);
			break;
		case 5:
			phiXY_Vx = MMCal.phi_sqV_5(sqVx, deltaXY);      	 
			wXY_Vx = MMCal.w_V_5(sqVx, deltaXY);
			break;
		case 6:
			phiXY_Vx = MMCal.phi_sqV_6(sqVx, deltaXY);      	 
			wXY_Vx = MMCal.w_V_6(sqVx, deltaXY);  
			break;
		case 7:
			phiXY_Vx = MMCal.phi_sqV_7(sqVx, deltaXY);      	 
			wXY_Vx = MMCal.w_V_7(sqVx, deltaXY);  
			break;
		}
		switch (phiZ) {
		case 1:
			phiZ_Vx = MMCal.phi_sqV_1(sqVtx, deltaZ);   
			wZ_Vx = MMCal.w_V_1(sqVtx, deltaZ); 
			break;
		case 2:
			phiZ_Vx = MMCal.phi_sqV_2(sqVtx, deltaZ);      	 
			wZ_Vx = MMCal.w_V_2(sqVtx,  deltaZ);
			break;
		case 3:
			phiZ_Vx = MMCal.phi_sqV_3(sqVtx,  deltaZ);      	 
			wZ_Vx = MMCal.w_V_3(sqVtx,  deltaZ);
			break;
		case 4:
			phiZ_Vx = MMCal.phi_sqV_4(sqVtx,  deltaZ);      	 
			wZ_Vx = MMCal.w_V_4(sqVtx,  deltaZ);
			break;
		case 5:
			phiZ_Vx = MMCal.phi_sqV_5(sqVtx,  deltaZ);      	 
			wZ_Vx = MMCal.w_V_5(sqVtx,  deltaZ);
			break;    	
		}

		double[][][] dphiXY_Vvx =  MMCal.dotMult(MMCal.copy(Vvx), wXY_Vx);
		double[][][] dphiXY_Vhx =   MMCal.dotMult(MMCal.copy(Vhx), wXY_Vx);

		double[][][] dphiZ_Vtx =  MMCal.dotMult(MMCal.copy(Vtx), wZ_Vx);
		double[][][] Hxy =  MMCal.minus(MMCal.copy(Hx),y);
		double J = 1.0/2.0 * MMCal.sum(MMCal.square(MMCal.concatAll(Hxy)))[0];
		//IJ.showMessage( "About optimisme...","phiXY_Vx "+toString(phiXY_Vx));
		//              "info  dans init process   PSFPreparator"+ String.valueOf(resIN[2]) + " \n" 
		//              		+ String.valueOf(resIN[1]) + " \n"  
		//              		+ String.valueOf(resIN[0]) + " \n"  
		//               		+ String.valueOf(this.dataDimPSF[2])  +" \n"
		//              		+ String.valueOf(this.dataDimPSF[1])  +" \n"
		//              		+ String.valueOf(this.dataDimPSF[0])  +" \n"); 


		//System.out.print("phiXY_Vx "+toString(phiXY_Vx));
		double[][][] lphiVx = MMCal.add( MMCal.dotMult(phiXY_Vx, lambdaXY) , MMCal.dotMult(phiZ_Vx,lambdaZ));
		double F =  J +  MMCal.sum(MMCal.concatAll(lphiVx))[0] ;
		double[] tabF = new double[1];
		double[][][] dJ = this.apply_PSFadjvar3D(Hxy);  
		double[][][] dphiXY_Vvx_op = MMCal.Vvoperatoradj(dphiXY_Vvx);   // dphiXY_Vvx corrupted
		dphiXY_Vvx = null;
		double[][][] dphiXY_Vhx_op = MMCal.Vhoperatoradj(dphiXY_Vhx); // dphiXY_Vhx  corrupted
		dphiXY_Vhx = null;
		double[][][] dphiXY_Vtx_op = MMCal.Vtoperatoradj(dphiZ_Vtx); // dphiZ_Vtx corrupted
		dphiZ_Vtx =null;
		double[][][] dF = new double[dphiXY_Vtx_op.length][dphiXY_Vtx_op[0].length][dphiXY_Vtx_op[0][0].length];
		for (int k = 0; k < dJ.length  ; k++)
			for (int j = 0 ; j < dJ[0].length ; j++)
				for ( int i = 0 ; i < dJ[0][0].length  ; i++)
				{
					dF[k][j][i] = dJ[k][j][i] + lambdaXY * ( dphiXY_Vvx_op[k][j][i] + dphiXY_Vhx_op[k][j][i])
							+ lambdaZ *  dphiXY_Vtx_op[k][j][i] ;
				}

		if(xmin != Integer.MIN_VALUE){
			int[][][] lte =  MMCal.lte(xmin, x);
			int[][][] lt =  MMCal.lt(xmin, x);
			double[] xred = MMCal.getWithFilter( x , lt);
			for (int i = 0 ; i < xred.length ; i++){
				xred[i] = 1.0/2.0* (xred[i] -  xmin); 
			}
			F = F + this.mmparam.eta * MMCal.norm(xred ,"fro");    	      	      	   
			dF = MMCal.add(dF,MMCal.dotMult(MMCal.dotMult(MMCal.add(MMCal.copy(x),-xmin) , lte), this.mmparam.eta));
		}
		if(xmax!= Integer.MAX_VALUE){ 		
			int[][][] gte = MMCal.gte(xmax, x);
			int[][][] gt = MMCal.gt(xmax, x);
			double[] xred =  MMCal.getWithFilter( x , gt);
			for (int i =0 ; i < xred.length ; i++)
				xred [i] = 1.0/2.0* (xred[i] - xmax);       		
			F = F + this.mmparam.eta * MMCal.norm(xred , "fro");
			dF = MMCal.add(dF, MMCal.dotMult(MMCal.dotMult(MMCal.add(MMCal.copy(x),-xmax) , gte), this.mmparam.eta));
		}
		tabF[0] = F;
		li.add(tabF);
		li.add(dF);
		li.add(wXY_Vx);
		li.add(wZ_Vx);
		return li;

	}


	/**
	 * @param A
	 * @return
	 */
	private static String toString(int[] A){
		String t ="[";
		for(int i=0 ; i < A.length ; i++){ 
			t = t +  String.valueOf(A[i]) + " ";
		}
		t = t + "]";
		return t;
	}
	protected static String toString(double[] A){
		String t ="[";
		for(int i=0 ; i < A.length ; i++){ 
			t = t +  String.valueOf(A[i]) + " ";
		}
		t = t + "]\n";
		return t;
	}
	protected static String toString(double[][] A){
		String t ="[";
		for(int j=0 ; j < A.length ; j++)
		{
			t = t + "[";
			for(int i=0 ; i < A[0].length ; i++){ 
				t = t +  String.valueOf(A[j][i]) + " ";
			}
			t = t + "]\n";
		}
		t = t + "]\n";
		return t;
	}
	protected static String toString(double[][][] A){
		String t ="[";
		for(int k=0 ; k < A.length ; k++){
			t = t + "[";	
			for(int j=0 ; j < A[0].length ; j++)
			{
				t = t + "[";
				for(int i=0 ; i < A[0][0].length ; i++){ 
					t = t +  String.valueOf(A[k][j][i]) + " ";
				}
				t = t + "]\n";
			}
			t = t + "]\n";
		}
		t = t + "]\n";
		return t;
	}	
	
	
	/**
	 * This Class implements attributs of algorithm parameters
	 * MMParameter encapsulates all different parameters of image and algorithms
     * Copyright:    Copyright (c) 2016
     * Company:      Labex Archimède
     *               AMU Aix Marseille Université
     * @author       Dominique Benielli
     * @author       Caraline Chaux-Moulin
     * @author       Sandrine Anthoine
     * @version 0.0 
	 *
	 */
	public static class MMParameter
	{
		public  int NbIt;   // visibility changed for test
		protected  double regul;
		protected  int T;
		protected  double regulZ;
		protected  int phixy;
		protected  int phiz; 
		protected  double TZ; 
		protected  int eta;
		protected  double[] resIM;
		protected double[] resPSF;
		protected String  unitRes;
		public MMParameter(int NbIt, 
				double regul, int T, double regulZ ,
				int phixy, int phiz, double TZ, 
				int eta, double[] resIM,double[] resPSF, String  unitRes ) 
		{ 
			this.NbIt = NbIt; 
			this.regul = regul;
			this.T = T;
			this.regulZ = regulZ;
			this.phixy = phixy;
			this.phiz = phiz;
			this.TZ = TZ;
			this.eta = eta;
			this.resIM = resIM;
			this.resPSF = resPSF;
			this.unitRes = unitRes;
		}
		
		public String toString(){

			return "NbIt : " + Integer.toString(this.NbIt) + " \n" +
					"regul : "+ Double.toString(this.regul) + " \n" +
					"T : "+ Integer.toString(this.T) + " \n" +
					"regulZ : " + Double.toString(this.regulZ) + " \n" +
					"phixy : " + Integer.toString(this.phixy) + " \n" +
					"phiz : " + Integer.toString(this.phiz) + " \n" +
					"TZ : " + Double.toString(this.TZ) + " \n" +
					"eta : " + Integer.toString(this.eta) + " \n" +
					"resIM : " + MM.toString(this.resIM) + " \n" +
					"resPSF : " + MM.toString(this.resPSF) + " \n" +
					"unitRes : "+ this.unitRes + " \n";

		}

	}
}
