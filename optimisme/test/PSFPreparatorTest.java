package optimisme.test;

import ij.ImagePlus;
import ij.ImageStack;
import ij.io.Opener;
import ij.process.ImageProcessor;

import optimisme.MM.MMParameter;
import optimisme.MMCal;
import optimisme.PSFPreparator;

import org.junit.After;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

public class PSFPreparatorTest {
	private ImagePlus impPsf = null;
	private ImagePlus impPSFCal  = null;
	private double[][][] dataPsfInit ;
	private double[][][] dataPsf;
	
	@Before
	public void setUpBefore() throws Exception {
	
		final String path =   System.getProperty("user.dir")+ "/data/";
		final String imageNamePsf = "psf_synth.tiff";
		final String imageNamePsfCal = "psf_synth_subsampled.tiff";
		Opener op = new Opener();
		this.impPsf = op.openImage(path, imageNamePsf);
		this.impPSFCal = op.openImage(path, imageNamePsfCal);
	}
	
	@After
	public  void tearDownAfter() throws Exception {
		this.impPsf.close();
		this.impPSFCal.close();
	//	this.imp.close();
	}	
	
	@Test
	public void fInterp3DTest(){	
    	int NbIt = 1000;
    	double regul = 1e-1;
    	int T = 1;
    	double regulZ = 1e-2;
    	double TZ  = 0.01;
    	int eta = 1;
    	int phixy = 4;
    	int phiz = 5;
    	int numberChannel = 5;
    	int numberChanneli = 1;
    	
    	double[][] expected ={{ 3.6884e-53  , 2.2774e-47 ,  2.3765e-42 ,
    		                    4.1916e-38  , 1.2495e-34 ,  6.2952e-32, 5.3606e-30 ,
    		                    7.7149e-29 ,  1.8766e-28 ,  7.7149e-29  , 5.3606e-30 ,  
    		                    6.2952e-32  ,1.2495e-34 ,  4.1916e-38 ,  2.3765e-42  , 
    		                    2.2774e-47},
    		                    { 2.2774e-47 ,  1.4061e-41 ,  1.4674e-36 ,  2.5881e-32 , 
    		                    	7.7149e-29  , 3.8869e-26,
    		                      3.3098e-24  , 4.7635e-23 ,  1.1587e-22  , 4.7635e-23  ,
    		                      3.3098e-24  , 3.8869e-26,
    		                    7.7149e-29 ,  2.5881e-32  , 1.4674e-36 ,  1.4061e-41},
    		                    {2.3765e-42 ,  1.4674e-36 ,  1.5313e-31  , 2.7008e-27  ,
    		                    	8.0509e-24,   4.0562e-21,
    		                    	3.4540e-19 ,  4.9709e-18  , 1.2091e-17  , 4.9709e-18 ,
    		                    	3.4540e-19 ,  4.0562e-21,
    		                    	 8.0509e-24 , 2.7008e-27 ,  1.5313e-31  , 1.4674e-36},
    		                    { 4.1916e-38  , 2.5881e-32 ,  2.7008e-27 , 4.7635e-23 , 
    		                    		 1.4200e-19 ,  7.1541e-17,
    		                    6.0919e-15 ,  8.7674e-14  , 2.1326e-13 ,  8.7674e-14  ,
    		                    6.0919e-15 ,  7.1541e-17,
    		                    1.4200e-19  , 4.7635e-23   ,2.7008e-27 ,  2.5881e-32	 
    		                    },
    		                    { 1.2495e-34 ,  7.7149e-29 ,  8.0509e-24  , 1.4200e-19  , 
    		                    	4.2329e-16 ,  2.1326e-13,
    		                    1.8160e-11 ,  2.6135e-10 ,  6.3572e-10 ,  2.6135e-10  ,
    		                    1.8160e-11 ,  2.1326e-13,
    		                    4.2329e-16 ,  1.4200e-19 ,  8.0509e-24  , 7.7149e-29
    		                    },
    		                    { 6.2952e-32  , 3.8869e-26 ,  4.0562e-21 ,  7.1541e-17 ,
    		                    2.1326e-13  , 1.0745e-10,
    		                    9.1492e-09 ,  1.3168e-07 ,  3.2029e-07 ,  1.3168e-07 ,
    		                    9.1492e-09 ,  1.0745e-10,
    		                    2.1326e-13 ,  7.1541e-17 ,  4.0562e-21 ,  3.8869e-26},
    		                    { 5.3606e-30 ,  3.3098e-24  , 3.4540e-19 ,  6.0919e-15 ,
    		                    	1.8160e-11 , 9.1492e-09,
    		                    7.7908e-07 ,  1.1212e-05 ,  2.7273e-05  , 1.1212e-05 ,
    		                    7.7908e-07 ,  9.1492e-09,
    		                    1.8160e-11  , 6.0919e-15  , 3.4540e-19 ,  3.3098e-24},
    		                    { 7.7149e-29 ,  4.7635e-23 ,  4.9709e-18 ,  8.7674e-14 , 
    		                    	2.6135e-10 ,  1.3168e-07,
    		                    1.1212e-05  , 1.6137e-04  , 3.9252e-04 ,  1.6137e-04  ,
    		                    1.1212e-05  , 1.3168e-07,
    		                     2.6135e-10 ,  8.7674e-14 ,  4.9709e-18 ,  4.7635e-23},
    		                     { 1.8766e-28 ,  1.1587e-22 ,  1.2091e-17  , 2.1326e-13 ,
    		                    6.3572e-10 ,  3.2029e-07,
    		                    2.7273e-05 ,  3.9252e-04 ,  9.5477e-04  , 3.9252e-04  ,
    		                    2.7273e-05,   3.2029e-07,
    		                   	6.3572e-10 ,  2.1326e-13 ,  1.2091e-17  , 1.1587e-22  },
    		                   	 {7.7149e-29 ,  4.7635e-23 ,  4.9709e-18 ,  8.7674e-14 ,
    		                   		2.6135e-10  , 1.3168e-07,
    		                   	 1.1212e-05  , 1.6137e-04  , 3.9252e-04 ,  1.6137e-04  ,
    		                   	 1.1212e-05 ,  1.3168e-07,
    		              	   2.6135e-10 ,  8.7674e-14 ,  4.9709e-18  , 4.7635e-23	},
    		              	   {  5.3606e-30,   3.3098e-24 ,  3.4540e-19  , 6.0919e-15 ,
    		              		   1.8160e-11 ,  9.1492e-09,
    		              		 7.7908e-07 ,  1.1212e-05 ,  2.7273e-05 ,  1.1212e-05 ,
    		              		 7.7908e-07 ,  9.1492e-09,
    		              	   1.8160e-11  , 6.0919e-15 ,  3.4540e-19 ,  3.3098e-24 },
    		              	   {6.2952e-32 ,  3.8869e-26 ,  4.0562e-21 ,  7.1541e-17  ,
    		              		   2.1326e-13  , 1.0745e-10,
    		              		 9.1492e-09 ,  1.3168e-07 ,  3.2029e-07 ,  1.3168e-07 ,
    		              		 9.1492e-09 ,  1.0745e-10,
    		              		2.1326e-13 ,  7.1541e-17 ,  4.0562e-21 ,  3.8869e-26 },
    		              		{1.2495e-34 ,  7.7149e-29 ,  8.0509e-24 ,  1.4200e-19  ,
    		              			4.2329e-16  , 2.1326e-13,
    		              		 1.8160e-11  , 2.6135e-10 ,  6.3572e-10 ,  2.6135e-10  , 
    		              		 1.8160e-11  , 2.1326e-13,
    		              		 4.2329e-16 ,  1.4200e-19 ,  8.0509e-24 ,  7.7149e-29 },
    		              		{4.1916e-38 ,  2.5881e-32 ,  2.7008e-27 ,  4.7635e-23 ,
    		              		1.4200e-19,  7.1541e-17,
    		              		6.0919e-15 ,  8.7674e-14 ,  2.1326e-13  , 8.7674e-14 , 
    		              		6.0919e-15  , 7.1541e-17,
    		              		1.4200e-19 ,  4.7635e-23  , 2.7008e-27  , 2.5881e-32 },
    		              		 {2.3765e-42 ,  1.4674e-36 ,  1.5313e-31 ,  2.7008e-27 , 
    		              			8.0509e-24 ,  4.0562e-21,
    		              		3.4540e-19 ,  4.9709e-18 ,  1.2091e-17 ,  4.9709e-18 ,
    		              		3.4540e-19 ,  4.0562e-21,
    		              		 8.0509e-24  , 2.7008e-27 ,  1.5313e-31 ,  1.4674e-36 },
    		              		{ 2.2774e-47,   1.4061e-41 ,  1.4674e-36  , 2.5881e-32  ,
    		              			 7.7149e-29 ,  3.8869e-26,
    		              		 3.3098e-24 ,  4.7635e-23  , 1.1587e-22 ,  4.7635e-23 , 
    		              		 3.3098e-24  , 3.8869e-26,
    		              		7.7149e-29  , 2.5881e-32 ,  1.4674e-36  , 1.4061e-41 } };
    	  	
    	double[] resPSF ={0.02, 0.02, 0.1};
    	double[] resIM ={0.08, 0.08, 0.2};
    	int[] sizeOUT = {Integer.MAX_VALUE, Integer.MAX_VALUE, Integer.MAX_VALUE};
    	this.dataPsfInit =  populateData(this.impPsf.getImageStack()); 
    	
    	double[][][] Hi = MMCal.fInterp3D(this.dataPsfInit, resPSF, resIM);
        double[][][] maximum = MMCal.max(Hi);
       // System.out.print("maximum :"+toString(maximum));
        for (int j =0 ; j< Hi[7].length; j++){
        	 for (int i =0 ; i< Hi[7][j].length; i++){
        	Assert.assertEquals(expected[j][i],Hi[7][j][i],10e-9);
        	 }
        }    
	}
	
	@Test
	public void processTest(){	
    	int NbIt = 1000;
    	double regul = 1e-1;
    	int T = 1;
    	double regulZ = 1e-2;
    	double TZ  = 0.01;
    	int eta = 1;
    	int phixy = 4;
    	int phiz = 5;
    	int numberChannel = 5;
    	int numberChanneli = 1;
    	double[] resPSF ={0.02, 0.02, 0.1};
    	double[] resIM ={0.08, 0.08, 0.2};
    	double[][][] expected;
    	int[] sizeOUT = {Integer.MAX_VALUE, Integer.MAX_VALUE, Integer.MAX_VALUE};
    	this.dataPsfInit =  populateData(this.impPsf.getImageStack()); 
    	this.dataPsf =  populateData(this.impPSFCal.getImageStack()); 
		MMParameter param = new MMParameter(NbIt, 
        		regul, T, regulZ ,
        		phixy, phiz, TZ, 
        		eta, resIM, resPSF, this.impPsf.getCalibration().getUnits() ) ;
	    PSFPreparator psfPreparator = new PSFPreparator(); 
	   
	    expected = psfPreparator.process(dataPsfInit, param, sizeOUT);

	    //    System.out.print(String.valueOf(this.dataPsf.length)+" "+String.valueOf(expected.length));
	//    System.out.print(String.valueOf(this.dataPsf[0].length)+" "+String.valueOf(expected[0].length));
	//    System.out.print(String.valueOf(this.dataPsf[0][0].length)+" "+String.valueOf(expected[0][0].length));
	    for (int k =0 ; k< this.dataPsf.length; k++){
        for (int j =0 ; j< this.dataPsf[0].length; j++){
       	 for (int i =0 ; i< this.dataPsf[0][0].length; i++){
       	 Assert.assertEquals(this.dataPsf[k][j][i], expected[k][j][i],10e-6);
       	 }
       }}
	}
	
    private double[][][] populateData(ImageStack source){
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

    
	private static String toString(double[][][] A){
		String t ="[";
		for(int k=0 ; k < A.length ; k++){
			t = t + "[";	
			for(int j=0 ; j < A[0].length ; j++)
			{
				t = t + "[";
				for(int i=0 ; i < A[0][0].length ; i++){ 
					t = t +  String.valueOf(A[k][j][i]) + " ";
				}
				t = t + "]";
		}
			t = t + "]";
		}
		t = t + "]";
		return t;
	}
	
}
