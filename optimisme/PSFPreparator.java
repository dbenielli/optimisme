package optimisme;

import ij.*;
import optimisme.MM.MMParameter;
import optimisme.MMCal;
import optimisme.MMCal.Squeeze;

/**
 * Title:        OPTIMIST --- A 3D  implementation for ImageJ
 * Copyright:    Copyright (c) 2017
 * Company:      Labex Archimède
 *               AMU Aix Marseille Université
 * @author       Dominique Benielli
 * @author       Caraline Chaux-Moulin
 * @author       Sandrine Anthoine
 * @version 1.2 
 */

public class PSFPreparator 
{
    private MMParameter mmparam;
    private int[] dataDim = new int[3];
    private int[] dataDimPSF = new int[3];

    public PSFPreparator ()
    {
      ;
    }

 
    /**
     * This Method is performed in the case of of the 
     * image input resolution is different as the PSF one.
     * The process method compute a new PSF with the same resolution of input image
     * @param psfdata This is 3D psf data
     * @param mmparam This is the MM parameters incuding images resolutions
     * @param sizeOut Table containing the image input dimension
     * @return This return the 3D  calculated PSF
     */
    public double[][][] process(double[][][] psfdata , MMParameter mmparam, int[] sizeOut  ) { 
    	this.dataDimPSF[2] = psfdata.length;
    	this.dataDimPSF[1] = psfdata[ 0 ].length;
    	this.dataDimPSF[0] = psfdata[ 0 ][ 0 ].length;

    	double[] resIN = mmparam.resPSF;
    	double[] resOUT = mmparam.resIM; 	
  //     	IJ.showMessage( "About optimisme...",
  //              "info  dans init process   PSFPreparator"+ String.valueOf(resIN[2]) + " \n" 
  //              		+ String.valueOf(resIN[1]) + " \n"  
  //              		+ String.valueOf(resIN[0]) + " \n"  
 //               		+ String.valueOf(this.dataDimPSF[2])  +" \n"
  //              		+ String.valueOf(this.dataDimPSF[1])  +" \n"
  //              		+ String.valueOf(this.dataDimPSF[0])  +" \n"); 
    	
    	IJ.showProgress(0.06);
    	double[][][] Hi = MMCal.fInterp3D(psfdata, resIN, resOUT);
    	
    	IJ.showProgress(0.3);
    	double[][][] Hib = MMCal.copy(Hi);   
    	IJ.showProgress(0.45);
    	double maxhi =  MMCal.max(MMCal.concatAll(Hib))[0];
   
    	Hib = MMCal.setValueWithFilter(0.0, Hib, MMCal.lt(0.5*maxhi,Hib));  	   	
    	Squeeze sq = new MMCal.Squeeze();
    	 	
    	
    	// ************************
    	// y0 calculation
    	// get the 2d matrix for y0
    	// *************************
        sq.cal(MMCal.sum(Hib,1));
        double[][] inter2 = (double[][]) sq.getSqueeze().get(0);
        sq.cal(MMCal.sum(MMCal.nonConjugateTranspose(inter2)));

        double[] inter1 = (double[]) sq.getSqueeze().get(0);                    	
        int y0 = (((int[]) MMCal.maxWithIndex(inter1).get(1)))[0];
        
 
    	// ************************
    	// x0 calculation       
    	// get the 2d matrix for x0
    	// *************************
        sq.cal(MMCal.sum(Hib,2)); 
        inter2 = (double[][]) sq.getSqueeze().get(0);
        sq.cal(MMCal.sum(inter2));
        inter1 = (double[]) sq.getSqueeze().get(0);        
        int x0 = (((int[]) MMCal.maxWithIndex(inter1).get(1)))[0];
 
        
    	// ************************
    	// z0 calculation             
    	// get the 2d matrix for z0
    	// *************************     
        sq.cal(MMCal.sum(Hib,0));
        inter2 = (double[][]) sq.getSqueeze().get(0);
        sq.cal(MMCal.sum(inter2));
        inter1 = (double[]) sq.getSqueeze().get(0);        
        int z0 = (((int[]) MMCal.maxWithIndex(inter1).get(1)))[0];       
        
        inter1 = null;
        inter2= null;
        
        //Matrix size
        int M2 = Hi[0].length;
        int M1 = Hi[0][0].length;
        int M3 = Hi.length;        
        
        //Center of 3D volume aroud PSD
        int rX = Math.min(Math.floorDiv(M1,2),(sizeOut[0]-1)/2);
        int rY = Math.min(Math.floorDiv(M2,2),(sizeOut[1]-1)/2);  
        int rZ = Math.min(Math.floorDiv(M3,2),(sizeOut[2]-1)/2); 
        
        int[] r0X = new int[x0+rX-(x0-rX)+1];
        for (int i=0 ; i< r0X.length ; i++)
        	r0X[i] = x0-rX + i;
        
        int[] r0Y = new int[y0+rY-(y0-rY)+1];
        for (int i=0 ; i< r0Y.length ; i++)
        	r0Y[i] = y0-rY + i;
        
        int[] r0Z = new int[z0+rZ-(z0-rZ)+1];
        for (int i=0 ; i< r0Z.length ; i++)
        	r0Z[i] = z0-rZ + i;
    
        int[] indX = MMCal.max(MMCal.min(r0X,M1-1),0);       
        int[] indY = MMCal.max(MMCal.min(r0Y,M2-1),0);
        int[] indZ = MMCal.max(MMCal.min(r0Z,M3-1),0);
     
        double[][][] H = MMCal.getIndex(Hi, indX, indY, indZ);

        // Remove negative values and normalize
        H = MMCal.max(H, 0.0);
        double normH = MMCal.sum( MMCal.concatAll(H))[0];
        H = MMCal.rDivide ( H, normH ); 
        
        //New Matrix size
      //  M2 = H[0].length;
      //  M1 = H[0][0].length;
      //  M3 = H.length;      
             
        return H;
        
        }

}
