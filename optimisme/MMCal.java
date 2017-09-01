package optimisme;
import ij.*;

import java.util.ArrayList;
import java.util.Arrays;
import edu.emory.mathcs.jtransforms.fft.DoubleFFT_3D;
import org.apache.commons.math3.analysis.interpolation.TricubicInterpolator;
import org.apache.commons.math3.analysis.interpolation.TricubicInterpolatingFunction;

/**
 * This Class provides some calculation methods to the MM algorithm
 * All of methods are based on the Matlab usual return
 * Copyright:    Copyright (c) 2017
 * Company:      Labex Archimède, I2M
 *               AMU Aix Marseille Université
 * @author       Dominique Benielli
 * @author       Caroline Chaux-Moulin
 * @author       Sandrine Anthoine
 * @version 1.2
 */
public  class MMCal{ 

	private static double maxValue = Double.NEGATIVE_INFINITY;

	/**
	 * Constructor
	 */
	public MMCal(){
	}

	/**
	 * 
	 * This stic method provide a deep copy of the input 3D data
	 * @param A double 3D element to copy
	 * @return The return is the deep copy of the input 3D data
	 */
	public static double[][][] copy(double[][][] A){
		double[][][] copy;
		copy = new double[A.length][][];
		for (int i = 0; i < A.length; i++) {
			copy[i] = new double[A[i].length][];
			for (int j = 0; j < A[i].length; j++) {
				copy[i][j] = new double[A[i][j].length];
				System.arraycopy(A[i][j], 0, copy[i][j], 0, 
						A[i][j].length);
			}
		}
		return copy;
	}

	/**
	 * This static method provide a deep copy of the input 2D data
	 * @param A double 2D element to copy
	 * @return The return is the deep copy of the input 2D data
	 */
	public static double[][] copy(double[][] A){
		double[][] copy;
		copy = new double[A.length][];
		for (int i = 0; i < A.length; i++) {
			copy[i] = new double[A[i].length];
			System.arraycopy(A[i], 0, copy[i], 0, 
					A[i].length);

		}
		return copy;
	}


	/**
	 * This static method realizes the setting of the value 'val' according 
	 * to the 'filter' containing 0 or 1 values, if filter value is 1 then
	 * the output data is setting with the value val, otherwise is inchanged
	 * @param val the output data setting to val where corresponding filter is equal to 1
	 * @param A input 3D data
	 * @param filter 3D data filled with 0 or 1 values
	 * @return 3D output data
	 */
	public static double[][][] setValueWithFilter(double val, 
			double[][][] A, 
			int[][][] filter ){
		double[][][] Hi = MMCal.copy(A);
		for (int k=0 ; k< Hi[0][0].length ; k++)
			for(int j=0 ; j< Hi[0].length ; j++)
				for(int i=0 ; i< Hi.length ; i++){
					if (filter[i][j][k] == 1)
						Hi[i][j][k] = val;
				}
		return Hi;
	}

	/**
	 * This static method return a 1D data of the input data 
	 * @param A The input 3D data
	 * @param filter The 3D filter setting with 0 or 1
	 * @return 1D output data filtered with filter
	 */
	public static double[] getWithFilter(double[][][] A, int[][][] filter ){
		double[] H = new double[A.length*A[0].length*A[0][0].length];
		double[] Hred = {};
		int index =0;
		for (int k=0 ; k< A.length ; k++)
			for(int j=0 ; j< A[0].length ; j++)
				for(int i=0 ; i< A[0][0].length ; i++){
					if (filter[k][j][i] == 1)
					{
						H[index]= A[k][j][i];
						index ++;
					}
				}
		if (index != 0) {
			Hred = Arrays.copyOf(H, index);
		}

		return Hred;
	}

	/**
	 * This static method set the value 'val' to the input 2D data 'V'
	 * @param V The input 2D data
	 * @param val The value to set
	 * @return The 2D output data
	 */
	public static  double[][] setValue(double[][] V, double val){
		for (int ai = 0; ai< V.length; ai++){
			Arrays.fill(V[ai], val);
		}
		return V;
	}

	/**
	 * This static method set the value 'val' to the input 3D data 'V'
	 * @param V The input 3D data
	 * @param val The value to set
	 * @return The 3D output data
	 */
	public static  double[][][] setValue(double[][][] V, double val){
		for (int ai = 0; ai < V.length; ai++)			
			for (int bi = 0; bi< V[0].length; bi++)
				Arrays.fill(V[ai][bi], val);			
		return V;
	}

	/**
	 * This static method set the value 'val' to the input 1D data 'V'
	 * @param V The input 1D data
	 * @param val The value to set
	 * @return The 1D output data
	 */
	public static  double[] setValue(double[] V, double val){
		//for (int bi = 0; bi < V.length; bi++)					
		//	V[bi] = val;
		Arrays.fill(V, val);
		return V;
	}

	/**
	 * The max static method
	 * @param A The 3D input data
	 * @param s The string 's' ("[]"), to indicate it is axis dependent
	 * @param a Integer to fix the axis of the maximum
	 * @return The 3D output data one of dimention is equal to 1
	 */
	public static double[][][] max(double[][][] A,String s , int a) 
	{
		return  (double[][][]) maxWithIndex(A, a).get(0);
	}

	/**
	 * The max static method
	 * @param A The 2D input data
	 * @param s The string 's' ("[]"), to indicate it is axis dependent
	 * @param a Integer to fix the axis of the maximum
	 * @return The 2D output data one of dimention is equal to 1
	 */public static double[][] max(double[][] A,String s , int a) 
	 {
		 return  (double[][]) maxWithIndex(A, a).get(0);
	 }

	 /**
	  * The max static method
	  * @param A The 1D input data
	  * @return The 1D output data the dimention is equal to 1
	  */
	 public static double[] max(double[] A){
		 double[] max = new double[1]; 
		 max[0] = maxValue;
		 for(int i=0 ; i< A.length ; i++){
			 if (A[i] > max[0])
				 max[0] = A[i];
		 }
		 return max;
	 }

	 /**
	  * The max static method, by default the axis is 0 if rows numbers is > 1
	  * @param A The 3D input data
	  * @return The 3D output data one of dimention is equal to 1
	  */public static double[][][] max(double[][][] A){
		  if (A[0].length > 1) return max(A,"[]", 0);
		  else 
			  if (A[0][0].length > 1) return max(A,"[]", 2);
			  else if  (A.length > 1) return max(A,"[]", 1);
			  else return A;
	  } 

	  /**
	   * This static method returns the input data wich are at least equal to 'val'
	   * @param A The 3D input data
	   * @param val The setting value at least
	   * @return The 3D data setted at least to val
	   */
	  public static double[][][] max(double[][][] A, double val){  	
		  double[][][] max = MMCal.copy(A);
		  for(int i=0 ; i< A.length ; i++)
			  for(int j=0 ; j< A[0].length ; j++)
				  for (int k=0 ; k< A[0][0].length ; k++)
					  if (A[i][j][k] < val){
						  max[i][j][k] = val;
					  }
		  return max;				
	  }

	  /**
	   * This static method concatenates all the value of input 
	   * columns are concatenated firt, after rows and slice
	   * @param A The 3D input data
	   * @return The 1D output data
	   */
	  public static double[] concatAll(double[][][] A){
		  double[] concat = new double[A.length*A[0].length*A[0][0].length];
		  int index1 = 0;
		  int index2 = 0;
		  for(int k=0 ; k< A.length ; k++){
			  for (int i=0 ; i< A[0][0].length ; i++){


				  for(int j=0 ; j< A[0].length ; j++)
					  concat[j+i*A[0].length+k*(A[0].length*A[0][0].length)] = A[k][j][i];
				  index1++;
			  }
			  index2++;
		  }
		  return concat;	
	  }

	  /**
	   * This static method return two elements, 
	   * the maximum values and the index of maximum values
	   * if a = 0
	   *   the maximum is performed along the columns axis 
	   * if a= 1
	   *   the maximum is performed along the slices axis 
	   * if a= 2
	   *   the maximum is performed along the rows axis      *   
	   * @param A the 3D input data
	   * @param a the axis to perform maximum (0,1,2)
	   * @return Two arrays , the maximum and the index
	   */
	  public static ArrayList maxWithIndex(double[][][] A, int a) {
		  ArrayList li = new ArrayList();
		  double[][][] max = new double[0][0][0];
		  int[][][] indexM = new int[0][0][0];

		  switch (a) {
		  case 1:
			  max = new double[1][A[0].length][A[0][0].length];
			  indexM = new int[1][A[0].length][A[0][0].length];
			  break;
		  case 0:
			  max = new double[A.length][1][A[0][0].length];
			  indexM = new int[A.length][1][A[0][0].length];
			  break;
		  case 2:
			  max = new double[A.length][A[0].length][1];  
			  indexM = new int[A.length][A[0].length][1]; 
			  break;
		  }
		  // INITIALIZATION VALUE OF maxValue
		  max = setValue(max, maxValue);

		  for(int i=0 ; i< A.length ; i++)
			  for(int j=0 ; j< A[0].length ; j++)
				  for (int k=0 ; k< A[0][0].length ; k++){
					  switch (a) {
					  case 1:
						  if (A[i][j][k] > max[0][j][k]){
							  max[0][j][k] = A[i][j][k];
							  indexM[0][j][k] = i;}
						  break;
					  case 0:
						  if (A[i][j][k] > max[i][0][k]){
							  max[i][0][k] = A[i][j][k];
							  indexM[i][0][k] = j;}                	
						  break;
					  case 2:
						  if (A[i][j][k] > max[i][j][0]){
							  max[i][j][0] = A[i][j][k];
							  indexM[i][j][0] = k;}
						  break;
					  }
				  }
		  li.add(0, max);
		  li.add(1, indexM);
		  return li;
	  } 

	  /**
	   * This static method return two elements, 
	   * the maximum values and the index of maximum values   *   
	   * @param A the 1D input data
	   * @return Two arrays of dimention 1 , the maximum and the index
	   */
	  public static ArrayList maxWithIndex(double[] A ) {
		  ArrayList li = new ArrayList();
		  double[] max = new double[1];
		  int[] indexM = new int[1];
		  // INITIALIZATION VALUE OF maxValue
		  max[0] = maxValue;
		  for(int i=0 ; i< A.length ; i++)
		  {
			  if (A[i] > max[0]){
				  max[0] = A[i];
				  indexM[0] = i;
			  }
		  }
		  li.add (0,max);
		  li.add(1,indexM);
		  return li;
	  }

	  /**
	   * This static method return two elements, 
	   * the maximum values and the index of maximum values
	   * if a = 0
	   *   the maximum is performed along the rows axis 
	   * if a= 1
	   *   the maximum is performed along the columns axis 
	   * @param A the 2D input data
	   * @param a the axis to perform maximum (0,1)
	   * @return Two arrays , the maximum and the index
	   */
	  public static ArrayList maxWithIndex(double[][] A, int a) 
	  {
		  ArrayList li = new ArrayList();
		  double[][] max = new double[0][0];
		  int[][] indexM = new int[0][0];
		  // INITIALIZATION of SHAPE
		  switch (a){
		  case 1:
			  max = new double[1][A[0].length];
			  indexM = new int[1][A[0].length];
			  break;
		  case 0:
			  max = new double[A.length][1];
			  indexM = new int[A.length][1];
			  break;
			  // INITIALIZATION VALUE OF maxValue	
		  }   
		  max = setValue(max, maxValue);
		  for(int i=0 ; i< A.length ; i++)
			  for(int j=0 ; j< A[0].length ; j++)
			  {
				  switch (a) {
				  case 1:
					  if (A[i][j] > max[0][j]){
						  max[0][j] = A[i][j];
						  indexM[0][j] = i;
					  }
					  break;
				  case 0:
					  if (A[i][j] > max[i][0]){
						  max[i][0] = A[i][j];
						  indexM[i][0] = j;
					  }
					  break;
				  }
			  }
		  li.add (0,max);
		  li.add(1,indexM);
		  return li;
	  }

	  /**
	   * This static method returns the input data wich are at least equal to 'val'
	   * @param A The 2D input data
	   * @param val The setting value at least
	   * @return The 2D data setted at least to val
	   */
	  public static double[][] max(double[][] A, double val){
		  double[][] max = MMCal.copy(A);
		  for(int i=0 ; i< A.length ; i++)
			  for(int j=0 ; j< A[0].length ; j++)
				  if (A[i][j] < val){
					  max[i][j] = val;
				  }
		  return max;
	  }

	  /**
	   * This static method returns the input data wich are at least equal to 'val'
	   * @param A The 1D input data
	   * @param val The setting value at least
	   * @return The 1D data setted at least to val
	   */
	  public static double[] max(double[] A, double val){	
		  double[] max =  Arrays.copyOf(A, A.length);
		  //new double[A.length];
		  //System.arraycopy(A, 0, max, 0, A.length);
		  for(int i=0 ; i< A.length ; i++)
			  if (A[i] < val){
				  max[i] = val;
			  }
		  return max;
	  }

	  /**
	   * This static method returns the input data wich are at least equal to 'val'
	   * @param A The 1D input data
	   * @param val The setting value at least
	   * @return The 1D data setted at least to val
	   */
	  public static int[] max(int[] A, int val){
		  int[] max =  Arrays.copyOf(A, A.length);
		  //	new int[A.length];
		  //System.arraycopy(A, 0, max, 0, A.length);
		  for(int i=0 ; i< A.length ; i++)
			  if (A[i] < val){
				  max[i] = val;
			  }
		  return max;
	  } 

	  /**
	   * This static method returns the input data wich are at most equal to 'val'
	   * @param A The 2D input data
	   * @param val The setting value at most
	   * @return The 2D data setted at most to val
	   */
	  public static double[] min(double[] A, double val){
		  double[] min = Arrays.copyOf(A, A.length);
		  for(int i=0 ; i< A.length ; i++)
			  if (A[i] > val){
				  min[i] = val;
			  }
		  return min;
	  }  

	  /**
	   * This static method returns the input data wich are at most equal to 'val'
	   * @param A The 1D input data
	   * @param val The setting value at most
	   * @return The 1D data setted at most to val
	   */
	  public static int[] min(int[] A, int val){
		  int[] min =  Arrays.copyOf(A, A.length);
		  for(int i=0 ; i< A.length ; i++)
			  if (A[i] > val){
				  min[i] = val;
			  }
		  return min;
	  }  

	  /**
	   * This static method convert input 2D double data to float 2D data
	   * @param datadouble2d input 2D double data
	   * @return 2D output float data
	   */
	  public static float[][] converttoFloat2d(double[][] datadouble2d){
		  float[][] data = new float[datadouble2d.length][datadouble2d[0].length];
		  for(int j =0 ; j< datadouble2d.length; j++ )
			  for(int i =0 ; i< datadouble2d[0].length; i++ )

				  data[j][i] = (float) datadouble2d[j][i];
		  return data;
	  }
	  
	  /**
	   * This static method convert input 2D double data to float 2D data
	   * and Reverse Dimension x and y for starck application
	   * @param datadouble2d input 2D double data
	   * @return 2D output float data
	   */
	  public static float[][] converttoFloat2dandReverseDim(double[][] datadouble2d){
		  float[][] data = new float[datadouble2d[0].length][datadouble2d.length];
		  for(int j =0 ; j< datadouble2d.length; j++ )
			  for(int i =0 ; i< datadouble2d[0].length; i++ )

				  data[i][j] = (float) datadouble2d[j][i];
		  return data;
	  }

	  /**
	   * This static method convert input 2D double data to float 1D data
	   * @param datadouble2d input 2D double data
	   * @return 1D output float data
	   */
	  public static float[] converttoFloat1d(double[][] datadouble2d){
		  float[] data = new float[datadouble2d.length * datadouble2d[0].length];
		  for(int j =0 ; j< datadouble2d.length; j++ )
			  for(int i =0 ; i< datadouble2d[0].length; i++ )

				  data[j*datadouble2d[0].length +i] = (float) datadouble2d[j][i];
		  return data;
	  }


	  /**
	   * The static method sum return the sum of the 3D input data
	   * by default the sum is along the columns (if row number of A > 1) 
	   * @param A The 3D input data
	   * @return the sum of A in one direction , the output 3D data has 1 dimension equal to one
	   */
	  public static double[][][] sum(double[][][] A){
		  if (A[0].length > 1) return sum(A, 0);
		  else 
			  if (A[0][0].length > 1) return sum(A, 1);
			  else if  (A.length > 1) return sum(A, 2);
			  else return A;
	  } 


	  /**
	   * The static method sum return the sum of the 3D input data
	   * if a = 0
	   *   the sum is performed along the columns axis 
	   * if a= 1
	   *   the sum is performed along the rows axis 
	   * if a= 2
	   *   the sum is performed along the slices axis    
	   * @param A The 3D input data
	   * @return the sum of A in one direction , the output 3D data has 1 dimension equal to one
	   */ 
	  public static double[][][] sum(double[][][] A, int a){
		  double[][][] sum = new double[0][0][0];
		  switch (a){
		  case 2:
			  sum = new double[1][A[0].length][A[0][0].length];
			  break;
		  case 0:
			  sum = new double[A.length][1][A[0][0].length];
			  break;
		  case 1:
			  sum = new double[A.length][A[0].length][1];
			  break;
		  }		
		  for(int k=0 ; k < A.length ; k++)
			  for(int j=0 ; j < A[0].length ; j++)
				  for (int i=0 ; i < A[0][0].length  ; i++){				
					  switch (a) {
					  case 2:
						  sum[0][j][i] += A[k][j][i];
						  break;
					  case 0:
						  sum[k][0][i] += A[k][j][i];
						  break;
					  case 1:
						  sum[k][j][0] += A[k][j][i];
						  break;
					  }
				  }
		  return sum;
	  }

	  /**
	   * The static method sum return the sum of the 2D input data
	   * by default the sum is along the columns (if row number of A > 1) 
	   * @param A The 2D input data
	   * @return the sum of A in one direction , the output 2D data has 1 dimension equal to one
	   */
	  public static double[][] sum(double[][] A){
		  if (A.length > 1) return sum(A, 0);
		  else 
			  if (A[0].length > 1) return sum(A, 1);
			  else return A;
	  }

	  /**
	   * The static method sum return the sum of the 2D input data
	   * if a = 0
	   *   the sum is performed along the columns axis 
	   * if a= 1
	   *   the sum is performed along the rows axis    
	   * @param A The 2D input data
	   * @return the sum of A in one direction , the output 2D data has 1 dimension equal to one
	   */ 
	  public static double[][] sum(double[][] A, int a){
		  double[][] sum = new double[0][0];
		  switch (a){
		  case 0:
			  sum = new double[1][A[0].length];
			  break;
		  case 1:
			  sum = new double[A.length][1];
			  break;
		  }		
		  for(int j=0 ; j < A.length ; j++)
			  for(int i=0 ; i < A[0].length ; i++){  				
				  switch (a) {
				  case 0:
					  sum[0][i] += A[j][i];
					  break;
				  case 1:
					  sum[j][0] += A[j][i];
					  break;
				  }

			  }
		  return sum;    	
	  }

	  /**
	   * The static method sum return the sum of the 1D input data  
	   * @param A The 1D input data
	   * @return output 1D data dimension equal to one
	   */ 
	  public static double[] sum(double[] A){
		  double[] sum = new double[1];
		  for(int i=0 ; i < A.length ; i++){ 
			  sum[0] += A[i];
		  }
		  return sum;	
	  }

	  /**
	   * This static method interCubic3 realizes the
	   * cubic 3D interpolation on the grid  (xi, yi, zi) of V on (x,y,z) grid
	   * @param z grid of the input data in slices direction
	   * @param y grid of the input data in columns direction
	   * @param x grid of the input data in rows direction
	   * @param zi grid of the output data in slices direction
	   * @param yi grid of the output data in columns direction
	   * @param xi grid of the output data in rows direction
	   * @param V input 3D data
	   * @return The 3D output data interpolated on grid (xi, yi, zi)
	   */
	  public static double[][][] interpCubic3(double[] z,double[] y,double[] x,
			  double[] zi,double[] yi,double[] xi,
			  double[][][] V ){

		  double[][][] Vi = new double[zi.length][yi.length][xi.length];   	
		  TricubicInterpolator tri = new TricubicInterpolator();			
		  TricubicInterpolatingFunction aa = tri.interpolate(z, y, x, V);

		  for (int k=0 ; k< zi.length ; k++)
			  for(int j=0 ; j< yi.length ; j++)
				  for(int i=0 ; i< xi.length ; i++)					 
					  Vi[k][j][i] = aa.value(zi[k], yi[j], xi[i]);

		  return Vi;
	  }

	  /**
	   * This static method fInterp3D realizes the
	   * cubic 3D interpolation in the new resolution
	   * @param V input 3D data
	   * @param resIN 3D tab of the input resolution (Input PSF resolution)
	   * @param resOUT 3D tab of the out resolution (output PSF resolution)
	   * @return The 3D output data interpolated
	   */
	  public static double[][][] fInterp3D(double[][][] V, 
			  double[] resIN, double[] resOUT)
	  {
		  int[] N = new int[3];
		  N[0] = V[0][0].length;
		  N[1] = V[0].length;
		  N[2] = V.length;

		  double taillex = resIN[0]*V[0][0].length;
		  double tailley = resIN[1]*V[0].length;
		  double taillez = resIN[2]*V.length;
		  double[] x = new  double[V[0][0].length];
		  double[] y = new  double[V[0].length];
		  double[] z = new  double[V.length];
		  int taillexi = (int) (Math.floor(resIN[0]*(double)(V[0][0].length-1)/resOUT[0]))+1;
		  int tailleyi = (int) (Math.floor(resIN[1]*(double)(V[0].length-1)/resOUT[1]))+1;
		  int taillezi = (int) (Math.floor(resIN[2]*(double)(V.length-1)/resOUT[2]))+1;
		  double[] xi = new  double[taillexi];
		  double[] yi = new  double[tailleyi];
		  double[] zi = new  double[taillezi];

		  for (int i=0; i < V[0][0].length; i++ ){
			  x[i] = i* resIN[0] ;    			
		  }
		  for (int i=0; i < V[0].length; i++ ){
			  y[i] = i* resIN[1] ;
		  }
		  for (int i=0; i < V.length; i++ ){
			  z[i] = i* resIN[2] ;
		  }

		  for (int i=0; i < taillexi; i++ ){
			  xi[i] = i* resOUT[0] ;    			
		  }
		  for (int i=0; i < tailleyi; i++ ){
			  yi[i] = i* resOUT[1] ;    			
		  }
		  for (int i=0; i < taillezi; i++ ){
			  zi[i] = i* resOUT[2] ;    			
		  }

		  //IJ.showProgress(0.1);
		  double[][][] Vi = interpCubic3(z,y,x,zi,yi,xi,V);
		  //IJ.showProgress(0.35);

		  //*double[][][] X;
		  //double[][][] Y;
		  //double[][][] Z;
		  //ArrayList al = this.meshgrid(x,y,z);
		  //X = al.get(0);
		  //Y = al.get(1);
		  //Z = al.get(2);
		  //double[][][] Xi;
		  //double[][][] Yi;
		  //double[][][] Zi;
		  //ArrayList ali = this.meshgrid(xi,yi,zi);
		  //Xi = ali.get(0);
		  //Yi = ali.get(1);
		  //Zi = ali.get(2);*/
		  return Vi;

	  }


	  /**
	   * This static method lt realizes the filtering of input 3D data
	   * LESS THAN the value Val
	   * @param val the threshold value for LESS THAN
	   * @param Hi input 3D data
	   * @return output 3D filtering data (with 0 or 1 values)
	   */
	  public static int[][][] lt(double val, double[][][] Hi){
		  int[][][] index = new int[Hi.length][Hi[0].length][Hi[0][0].length];
		  for(int i=0 ; i< Hi.length ; i++)
			  for(int j=0 ; j< Hi[0].length ; j++)
				  for (int k=0 ; k< Hi[0][0].length ; k++)
				  {
					  if (Hi[i][j][k] < val)
						  index[i][j][k] = 1;
				  }
		  return index;
	  }

	  /**
	   * This static method lte realizes the filtering of input 3D data
	   * LESS THAN or EQUAL the value Val
	   * @param val the threshold value for LESS THAN or EQUAL 
	   * @param Hi input 3D data
	   * @return output 3D filtering data (with 0 or 1 values)
	   */
	  public static int[][][] lte(double val, double[][][] Hi){
		  int[][][] index = new int[Hi.length][Hi[0].length][Hi[0][0].length];
		  for(int i=0 ; i< Hi.length ; i++)
			  for(int j=0 ; j< Hi[0].length ; j++)
				  for (int k=0 ; k< Hi[0][0].length ; k++)
				  {
					  if (Hi[i][j][k] <= val)
						  index[i][j][k] = 1;
				  }
		  return index;
	  }

	  /**
	   * This static method gt realizes the filtering of input 3D data
	   * GREAT THAN the value Val
	   * @param val the threshold value for GREAT THAN 
	   * @param Hi input 3D data
	   * @return output 3D filtering data (with 0 or 1 values)
	   */
	  public static int[][][] gt(double val, double[][][] Hi){
		  int[][][] index = new int[Hi.length][Hi[0].length][Hi[0][0].length];
		  for(int i=0 ; i< Hi.length ; i++)
			  for(int j=0 ; j< Hi[0].length ; j++)
				  for (int k=0 ; k< Hi[0][0].length ; k++)
				  {
					  if (Hi[i][j][k] > val)
						  index[i][j][k] = 1;
				  }
		  return index;
	  }  

	  /**
	   * This static method gte realizes the filtering of input 3D data
	   * GREAT THAN or EQUAL the value Val
	   * @param val the threshold value for GREAT THAN or EQUAL
	   * @param Hi input 3D data
	   * @return output 3D filtering data (with 0 or 1 values)
	   */
	  public static int[][][] gte(double val, double[][][] Hi){
		  int[][][] index = new int[Hi.length][Hi[0].length][Hi[0][0].length];
		  for(int i=0 ; i< Hi.length ; i++)
			  for(int j=0 ; j< Hi[0].length ; j++)
				  for (int k=0 ; k< Hi[0][0].length ; k++)
				  {
					  if (Hi[i][j][k] >= val)
						  index[i][j][k] = 1;
				  }
		  return index;
	  }  


	  /**
	   * This static method nonConjugateTranspose realizes the transposition non conjugate
	   * of input 3D data a(i,j) -> a(j,i)
	   * @param A input 3D data
	   * @return  output 3D data transposed
	   */
	  public static double[][] nonConjugateTranspose(double[][] A){
		  double[][] trans = new double[A[0].length][A.length];
		  for(int j=0 ; j< A.length ; j++)
			  for(int i=0 ; i< A[0].length ; i++){
				  trans[i][j] = A[j][i];
			  }
		  return trans;
	  }

	  /**
	   * This static method getIndex returns 3D data on the specified index
	   * @param A
	   * @param indexX columns index needed
	   * @param indexY rows index needed
	   * @param indexZ slices index needed
	   * @return output 3D data on indexX, indexY, indexZ
	   */
	  public static double[][][] getIndex(double[][][] A, int[] indexX, int[] indexY,int[] indexZ){
		  double[][][] Aindex= new double[indexZ.length][indexY.length][indexX.length];
		  for (int k=0 ; k< indexZ.length; k++)
			  for(int j=0 ; j< indexY.length ; j++)
				  for(int i=0 ; i< indexX.length ; i++){
					  Aindex[k][j][i] = A[indexZ[k]][indexY[j]][indexX[i]];
				  }
		  return Aindex;

	  }


	  /**
	   * This static method rDivide realizes the division of 3D input data by the specified value 'val'
	   * @param A input 3D data
	   * @param value specified value for the division
	   * @return output 3D data
	   */
	  public static double[][][] rDivide(double[][][] A, double value){
		  for (int k=0 ; k< A.length; k++)
			  for(int j=0 ; j< A[0].length ; j++)
				  for(int i=0 ; i< A[0][0].length ; i++){
					  A[k][j][i] = A[k][j][i]/value;
				  }
		  return A;

	  }


	  /**
	   * DON'T USED FUNCTION realizes convolution in the case 
	   * of function_handle = true , to be implemented
	   * @param D 2D input data from image
	   * @param ps 2D PSF from image
	   * @return 2D output data 
	   */
	  public static double[][] conv2Dadjoint_fourier(double[][] D, double[][] ps){
		  //	int p = Math.floorDiv(arg0, arg1)
		  //TO DO 	;
		  /*p=floor(size(a)/2);

	[N1,N2]=size(D);

	Dext = zeros(size(D)+2*p);
	Dext(p(1)+(1:N1),p(2)+(1:N2)) = D;

	A = conj(freqfilt2D(a,2*p(1)+N1,2*p(2)+N2));

	bF2 = ifft2(fft2(Dext).*A);
	bF2 = bF2(p(1)+(1:N1),p(2)+(1:N2));
		   * 	
		   */
		  return D;
	  }

	  /**
	   * This static method conv3D_fourier realizes convolution in the fourier space of inputs data
	   * @param D 3D input data from image
	   * @param ps 3D PSF from image
	   * @return 3D output data 
	   */
	  public static double[][][] conv3D_fourier(double[][][] D, double[][][] ps){
		  int[] p = new int[3];
		  int N1 = D[0][0].length;
		  int N2 = D[0].length;
		  int N3 = D.length;
		  p[0] = Math.floorDiv(ps[0][0].length, 2);
		  p[1] = Math.floorDiv(ps[0].length, 2);
		  p[2] = Math.floorDiv(ps.length, 2);
		  double[][][] Dext = new double[N3+2*p[2]][N2+2*p[1]][N1+2*p[0]];
		  for (int k = 0; k < N3  ; k++)
			  for (int j = 0; j < N2 ; j++)
				  for (int i = 0 ; i < N1  ; i++)
					  Dext[k+p[2]][j+p[1]][i+p[0]] = D[k][j][i];
		  // is a double complex
		  //   System.out.print("2p+N1 : "+ String.valueOf(2*p[0]+N1));
		  //   System.out.print("2p+N2 : "+ String.valueOf(2*p[1]+N2));
		  //   System.out.print("2p+N3 : "+ String.valueOf(2*p[2]+N3));
		  double[][][] Acomplex = freqfilt3D(ps,2*p[0]+N1,2*p[1]+N2,2*p[2]+N3);
		  //    System.out.print("Aextcomplex :"+toString(Acomplex)); 
		  //   
		  DoubleFFT_3D FFTCal1 = new DoubleFFT_3D(Dext.length,Dext[0].length,Dext[0][0].length);
		  double[][][] Dextcomplex = complexCopy(Dext);
		  FFTCal1.complexForward(Dextcomplex);


		  Dextcomplex = MMCal.dotMultComplex(Dextcomplex, Acomplex);
		  //      System.out.print("Dextcomplex mult:"+toString(Dextcomplex));  
		  Acomplex = null;

		  DoubleFFT_3D FFTCal2 = new DoubleFFT_3D(Dext.length,Dext[0].length,Dext[0][0].length);
		  FFTCal2.complexInverse(Dextcomplex, true);

		  //   System.out.print("Dextcomplex ifft:"+toString(Dextcomplex)); 

		  int[] indX = new int[N1];
		  for (int i = 0 ; i < N1 ; i++)
			  indX[i] = p[0] + i;
		  int[] indY = new int[N2];
		  for (int i = 0 ; i < N2 ; i++)
			  indY[i] = p[1] + i;   	
		  int[] indZ = new int[N3];
		  for (int i = 0 ; i < N3 ; i++)
			  indZ[i] = p[2] + i;
		  double[][][] bF2double = toDoubleReal( Dextcomplex);


		  Dextcomplex = null;
		  bF2double = MMCal.getIndex(bF2double, indX, indY, indZ);
		  return bF2double;
	  }

	  /**
	   * This static method  conj returns the complex conjugate of input and Complex 3D data
	   * the input array must be of size slices by rows by 2*columns, 
	   * physical layout of the input data is as follows:
	   * <pre>
	   * a[k3][k2][2*k1] = Re[k3][k2][k1], a[k3][k2][2*k1+1] = Im[k3][k2][k1],
	   * 0&lt;=k3&lt;slices, 0&lt;=k2&lt;columns, 0&lt;=k1&lt;rows,
	   * </pre>
	   * @param aComplex 3D and Complex input data
	   * @return 3D and Complex output data conjugate
	   */
	  public static double[][][] conj(double[][][] aComplex){
		  for (int k = 0; k < aComplex.length  ; k++)
			  for (int j = 0 ; j < aComplex[0].length ; j++)
				  for (int i = 0 ; i < aComplex[0][0].length/2  ; i++)
					  aComplex[k][j][2*i+1]= -aComplex[k][j][2*i+1];
		  return aComplex;
	  }

	  /**
	   * This static method  freqfilt3D returns the complex frequency output 3D data
	   * the output array is of size slices by rows by 2*columns, 
	   * physical layout of the output data is as follows:
	   * <pre>
	   * a[k3][k2][2*k1] = Re[k3][k2][k1], a[k3][k2][2*k1+1] = Im[k3][k2][k1],
	   * 0&lt;=k3&lt;slices, 0&lt;=k2&lt;columns, 0&lt;=k1&lt;rows,
	   * </pre>
	   * @param a  3D input data
	   * @param n Integer for L1 rows dimension
	   * @param m Integer for L2 columns dimension
	   * @param p Integer for L3 slices dimension
	   * @return 3D and Complex output frequency data 
	   */
	  public static double[][][] freqfilt3D(double[][][] a, int n, int m, int p){
		  int L1 = a[0][0].length;
		  int L2 = a[0].length;
		  int L3 = a.length;
		  double[][][] T  = new double[p][m][n];
		  int L13 = (int) Math.ceil(L1/2.0);
		  int L23 = (int) Math.ceil(L2/2.0);
		  int L33 = (int) Math.ceil(L3/2.0);

		  for (int k = L33-1; k < L3  ; k++)
		  {
			  for (int j = L23-1 ; j < L2 ; j++)
			  {
				  for (int i = 0 ; i < L13-1  ; i++)    
					  T[k-L33+1][j-L23+1][i-L13+n+1] = a[k][j][i];
				  for (int i = L13-1 ; i < L1  ; i++)
					  T[k-L33+1][j-L23+1][i-L13+1] = a[k][j][i];
			  }
			  for (int j = 0 ; j < L23-1 ; j++)
			  {
				  for (int i = 0 ; i < L13-1  ; i++)
					  T[k-L33+1][j-L23+m+1][i-L13+n+1] = a[k][j][i]; 
				  for (int i = L13-1 ; i < L1  ; i++)
					  T[k-L33+1][j-L23+m+1][i-L13+1] = a[k][j][i];
			  }

		  }

		  for (int k = 0; k < L33-1  ; k++)
		  {
			  for (int j = L23-1 ; j < L2 ; j++)
			  {
				  for (int i = 0 ; i < L13-1  ; i++) 
					  T[k-L33+p+1][j-L23+1][i-L13+n+1] = a[k][j][i];
				  for (int i = L13-1 ; i < L1  ; i++)    
					  T[k-L33+p+1][j-L23+1][i-L13+1] = a[k][j][i];
			  }
			  for (int j = 0 ; j < L23-1 ; j++)
			  {
				  for (int i = 0 ; i < L13-1  ; i++)
					  T[k-L33+p+1][j-L23+m+1][i-L13+n+1] = a[k][j][i];
				  for (int i = L13-1 ; i < L1  ; i++)
					  T[k-L33+p+1][j-L23+m+1][i-L13+1] = a[k][j][i]; 
			  }
		  }

		  // ComplexNum[][][] Tcomplex = toComplexNum(MMCal.reverseCube(T));
		  //  ComplexNum[][][] Tcomplex = toComplexNum(T);

		  DoubleFFT_3D FFTCal = new DoubleFFT_3D(T.length,T[0].length,T[0][0].length);
		  double[][][] Tcomplex = complexCopy(T);
		  FFTCal.complexForward(Tcomplex);

		  //  DoublePrecFFT3D FFTCal = new DoublePrecFFT3D(Tcomplex);
		  //  FFTCal.fft();
		  //  Tcomplex = FFTCal.getData(); 
		  return 		Tcomplex;
	  }


	  /**
	   * This static method  complexCopy returns the complex copy input 3D data
	   * the output array is of size slices by rows by 2*columns, 
	   * physical layout of the output data is as follows:
	   * <pre>
	   * a[k3][k2][2*k1] = Re[k3][k2][k1], a[k3][k2][2*k1+1] = Im[k3][k2][k1],
	   * 0&lt;=k3&lt;slices, 0&lt;=k2&lt;columns, 0&lt;=k1&lt;rows,
	   * </pre>
	   * @param T 3D real input data
	   * @return 3D Complex data
	   */
	  public static double[][][] complexCopy(double[][][] T){
		  double[][][] Tcomplex = new double[T.length][T[0].length][T[0][0].length * 2];
		  for(int k = 0 ; k< T.length ; k++)
			  for(int j = 0 ; j< T[0].length ; j++)
				  for(int i = 0 ; i< T[0][0].length ; i++){   			  				
					  Tcomplex[k][j][2*i] = T[k][j][i];    				
					  Tcomplex[k][j][2*i+1] =  0;
				  }
		  //	System.arraycopy(T[k][j], 0, Tcomplex[k][j], 0, T[k][j].length);
		  return Tcomplex;
	  }

	  /**
	   * This static method toString returns a string represention of 1D input data
	   * @param A 1D input data
	   * @return A string representation of input data
	   */
	  public static String toString(double[] A){
		  String t ="[";
		  for(int i=0 ; i < A.length ; i++){ 
			  t = t +  String.valueOf(A[i]) + " ";
		  }
		  t = t + "]\n";
		  return t;
	  }

	  /**
	   * This static method toString returns a string represention of 3D input data
	   * @param A 3D input data
	   * @return A string representation of input data
	   */
	  public static String toString(double[][][] A){
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
		  t = t + "]";
		  return t;
	  }	  


	  /**
	   * This static method  toDoubleReal  returns the real part of input 3D data
	   * the input array must be of size slices by rows by 2*columns, 
	   * physical layout of the input data is as follows:
	   * <pre>
	   * a[k3][k2][2*k1] = Re[k3][k2][k1], a[k3][k2][2*k1+1] = Im[k3][k2][k1],
	   * 0&lt;=k3&lt;slices, 0&lt;=k2&lt;columns, 0&lt;=k1&lt;rows,
	   * </pre> 
	   * @param c input complex 3D data
	   * @return The real part
	   */
	  public static double[][][] toDoubleReal(double[][][] c){
		  double[][][] re = new double[c.length][c[0].length][c[0][0].length/2];
		  for (int k = 0; k < re.length  ; k++)
			  for (int j = 0 ; j < re[0].length ; j++)
				  for (int i = 0 ; i < re[0][0].length  ; i++)
					  re[k][j][i] =  c[k][j][2*i];
		  return re;     
	  }    

	  /**
	   * This static method  sqrt caculs the sqrt of the input 3D data
	   * @param x input 3D data
	   * @return the sgrt of data
	   */
	  public static double[][][] sqrt(double[][][] x){
		  for (int k = 0; k < x.length  ; k++)
			  for (int j = 0 ; j < x[0].length ; j++)
				  for ( int i = 0 ; i < x[0][0].length  ; i++)
					  x[k][j][i] = Math.sqrt(x[k][j][i] );
		  return x;
	  }

	  /**
	   * This static method  log  caculs the logaritm of the input 3D data
	   * @param x input 3D data
	   * @return the logarithm of data
	   */
	  public static double[][][] log(double[][][] x){
		  for (int k = 0; k < x.length  ; k++)
			  for (int j = 0 ; j < x[0].length ; j++)
				  for ( int i = 0 ; i < x[0][0].length  ; i++)
					  x[k][j][i] = Math.log(x[k][j][i] );
		  return x;
	  }

	  /**
	   * This static method   caculs the power of  <code>d</code> 
	   * of the input 3D data
	   * @param x input 3D data
	   * @param d value of the power of
	   * @return the power <code>d</code>  of data
	   */
	  public static double[][][] pow(double[][][] x, double d){
		  for (int k = 0; k < x.length  ; k++)
			  for (int j = 0 ; j < x[0].length ; j++)
				  for ( int i = 0 ; i < x[0][0].length  ; i++)
					  x[k][j][i] = Math.pow(x[k][j][i], d);
		  return x;
	  }

	  /**
	   * This static method caculs the square 
	   * of the input 3D data	  
	   * @param x input 3D data
	   * @return 3D output : the square of the input 3D data
	   */
	  public static double[][][] square(double[][][] x){
		  for (int k = 0; k < x.length  ; k++)
			  for (int j = 0 ; j < x[0].length ; j++)
				  for ( int i = 0 ; i < x[0][0].length  ; i++)
					  x[k][j][i] = Math.pow(x[k][j][i],2);
		  return x;
	  }

	  /**
	   * This static method caculs the square 
	   * of the input 1D data	  
	   * @param x input 1D data
	   * @return 1D output : the square of the input 1D data
	   */
	  public static double[] square(double[] x){
		  for (int k = 0; k < x.length  ; k++)
			  x[k] = Math.pow(x[k],2);
		  return x;
	  }

	  /**
	   * This static method corresponds to 
	   * the option 1  of the phi function
	   * the penalty function 
	   * @param sqVx square of Vx input 3D data
	   * @param deltaXY double parameter
	   * @return The phi function of Vx square
	   */
	  public static double[][][] phi_sqV_1(double[][][] sqVx, double deltaXY){
		  double[][][] re = new double[sqVx.length][sqVx[0].length][sqVx[0][0].length];
		  for (int k = 0; k < sqVx.length  ; k++)
			  for (int j = 0 ; j < sqVx[0].length ; j++)
				  for ( int i = 0 ; i < sqVx[0][0].length  ; i++)
					  re[k][j][i] = 1.0 - Math.exp(- sqVx[k][j][i] / (2.0*Math.pow(deltaXY,2)));
		  return re;
	  }

	  /**
	   * This static method corresponds to 
	   * the option 2 of the phi function
	   * the penalty function 
	   * @param sqVx square of Vx input 3D data
	   * @param deltaXY double parameter
	   * @return The phi function of Vx square
	   */
	  public static double[][][] phi_sqV_2(double[][][] sqVx, double deltaXY){
		  double[][][] re = new double[sqVx.length][sqVx[0].length][sqVx[0][0].length];
		  for (int k = 0; k < sqVx.length  ; k++)
			  for (int j = 0 ; j < sqVx[0].length ; j++)
				  for ( int i = 0 ; i < sqVx[0][0].length  ; i++)
					  re[k][j][i] = sqVx[k][j][i] / ( (2.0*Math.pow(deltaXY,2)) + sqVx[k][j][i]);
		  return re;
	  }  

	  /**
	   * This static method corresponds to 
	   * the option 3  of the phi function
	   * the penalty function 
	   * @param sqVx square of Vx input 3D data
	   * @param deltaXY double parameter
	   * @return The phi function of Vx square
	   */
	  public static double[][][] phi_sqV_3(double[][][] sqVx, double deltaXY){
		  double[][][] re = new double[sqVx.length][sqVx[0].length][sqVx[0][0].length];
		  for (int k = 0; k < sqVx.length  ; k++)
			  for (int j = 0 ; j < sqVx[0].length ; j++)
				  for ( int i = 0 ; i < sqVx[0][0].length  ; i++)
					  re[k][j][i] = Math.log(1.0 + sqVx[k][j][i] / Math.pow(deltaXY,2));

		  return re;
	  } 

	  /**
	   * This static method corresponds to 
	   * the option 4  of the phi function
	   * the penalty function 
	   * @param sqVx square of Vx input 3D data
	   * @param deltaXY double parameter
	   * @return The phi function of Vx square
	   */
	  public static double[][][] phi_sqV_4(double[][][] sqVx, double deltaXY){
		  double[][][] re = new double[sqVx.length][sqVx[0].length][sqVx[0][0].length];
		  for (int k = 0; k < sqVx.length  ; k++)
			  for (int j = 0 ; j < sqVx[0].length ; j++)
				  for ( int i = 0 ; i < sqVx[0][0].length  ; i++)
					  re[k][j][i] = Math.sqrt(1.0 + sqVx[k][j][i] / Math.pow(deltaXY,2))  -1.0;                  
		  return re;
	  } 

	  /**
	   * This static method corresponds to 
	   * the option 5  of the phi function
	   * the penalty function 
	   * @param sqVx square of Vx input 3D data
	   * @param deltaXY double parameter
	   * @return The phi function of Vx square
	   */
	  public static double[][][] phi_sqV_5(double[][][] sqVx, double deltaXY){
		  double[][][] re = new double[sqVx.length][sqVx[0].length][sqVx[0][0].length];
		  for (int k = 0; k < sqVx.length  ; k++)
			  for (int j = 0 ; j < sqVx[0].length ; j++)
				  for ( int i = 0 ; i < sqVx[0][0].length  ; i++)
					  re[k][j][i] = (1.0/2.0) * sqVx[k][j][i] ;             
		  return re;
	  } 

	  /**
	   * This static method corresponds to 
	   * the option 6  of the phi function
	   * the penalty function 
	   * @param sqVx square of Vx input 3D data
	   * @param deltaXY double parameter
	   * @return The phi function of Vx square
	   */
	  public static double[][][] phi_sqV_6(double[][][] sqVx, double deltaXY){
		  double[][][] re = new double[sqVx.length][sqVx[0].length][sqVx[0][0].length];
		  for (int k = 0; k < sqVx.length  ; k++)
			  for (int j = 0 ; j < sqVx[0].length ; j++)
				  for ( int i = 0 ; i < sqVx[0][0].length  ; i++)
					  re[k][j][i] = 1.0 - Math.exp( - ( Math.sqrt( 1.0 + sqVx[k][j][i] / Math.pow(deltaXY,2)) -1.0) );
		  return re;
	  }
	  /**
	   * This static method corresponds to 
	   * the option 7  of the phi function,
	   * the penalty function 
	   * @param sqVx square of Vx input 3D data
	   * @param deltaXY double parameter
	   * @return The phi function of Vx square
	   */
	  public static double[][][] phi_sqV_7(double[][][] sqVx, double deltaXY){
		  double[][][] re = new double[sqVx.length][sqVx[0].length][sqVx[0][0].length];
		  final double  pow = 0.25;
		  for (int k = 0; k < sqVx.length  ; k++)
			  for (int j = 0 ; j < sqVx[0].length ; j++)
				  for ( int i = 0 ; i < sqVx[0][0].length  ; i++)
					  re[k][j][i] = Math.pow((1.0 + sqVx[k][j][i] /Math.pow(deltaXY,2)), pow )-1.0;
		  return re;
	  }

	  /**
	   * This static method corresponds to 
	   * the option 1  of  of the w function,
	   * the penalty function 
	   * @param sqVx square of Vx input 3D data
	   * @param deltaXY double parameter
	   * @return The phi function of Vx square
	   */
	  public static double[][][] w_V_1(double[][][] sqVx, double deltaXY){
		  double[][][] re = new double[sqVx.length][sqVx[0].length][sqVx[0][0].length];
		  for (int k = 0; k < sqVx.length  ; k++)
			  for (int j = 0 ; j < sqVx[0].length ; j++)
				  for ( int i = 0 ; i < sqVx[0][0].length  ; i++)
					  re[k][j][i] = (1.0/(Math.pow(deltaXY,2))) * Math.exp(- sqVx[k][j][i] /(2.0*Math.pow(deltaXY,2))); 
		  return re;
	  }

	  /**
	   * This static method corresponds to 
	   * the option 2  of  of the w function,
	   * the penalty function 
	   * @param sqVx square of Vx input 3D data
	   * @param deltaXY double parameter
	   * @return The phi function of Vx square
	   */
	  public static double[][][] w_V_2(double[][][] sqVx, double deltaXY){
		  double[][][] re = new double[sqVx.length][sqVx[0].length][sqVx[0][0].length];
		  for (int k = 0; k < sqVx.length  ; k++)
			  for (int j = 0 ; j < sqVx[0].length ; j++)
				  for ( int i = 0 ; i < sqVx[0][0].length  ; i++)
					  re[k][j][i] = (4.0*Math.pow(deltaXY,2) ) / Math.pow((2.0*Math.pow(deltaXY,2))+ sqVx[k][j][i], 2 ); 
		  return re;
	  }

	  /**
	   * This static method corresponds to 
	   * the option 3  of  of the w function,
	   * the penalty function 
	   * @param sqVx square of Vx input 3D data
	   * @param deltaXY double parameter
	   * @return The phi function of Vx square
	   */
	  public static double[][][] w_V_3(double[][][] sqVx, double deltaXY){
		  double[][][] re = new double[sqVx.length][sqVx[0].length][sqVx[0][0].length];
		  for (int k = 0; k < sqVx.length  ; k++)
			  for (int j = 0 ; j < sqVx[0].length ; j++)
				  for ( int i = 0 ; i < sqVx[0][0].length  ; i++)
					  re[k][j][i] = 2.0 / ( Math.pow(deltaXY,2) + sqVx[k][j][i]);
		  return re;
	  }

	  /**
	   * This static method corresponds to 
	   * the option 4  of  of the w function,
	   * the penalty function 
	   * @param sqVx square of Vx input 3D data
	   * @param deltaXY double parameter
	   * @return The phi function of Vx square
	   */
	  public static double[][][] w_V_4(double[][][] sqVx, double deltaXY){
		  double[][][] re = new double[sqVx.length][sqVx[0].length][sqVx[0][0].length];
		  for (int k = 0; k < sqVx.length  ; k++)
			  for (int j = 0 ; j < sqVx[0].length ; j++)
				  for ( int i = 0 ; i < sqVx[0][0].length  ; i++)
					  re[k][j][i] = (1/Math.pow(deltaXY,2)) / Math.sqrt(1 + sqVx[k][j][i] / Math.pow(deltaXY,2));
		  return re;
	  }

	  /**
	   * This static method corresponds to 
	   * the option 5  of  of the w function,
	   * the penalty function 
	   * @param sqVx square of Vx input 3D data
	   * @param deltaXY double parameter
	   * @return The phi function of Vx square
	   */
	  public static double[][][] w_V_5(double[][][] sqVx, double deltaXY){
		  double[][][] re = new double[sqVx.length][sqVx[0].length][sqVx[0][0].length];
		  MMCal.setValue(re, 1);
		  return re;
	  }

	  /**
	   * This static method corresponds to 
	   * the option 6  of  of the w function,
	   * the penalty function 
	   * @param sqVx square of Vx input 3D data
	   * @param deltaXY double parameter
	   * @return The phi function of Vx square
	   */
	  public static double[][][] w_V_6(double[][][] sqVx, double deltaXY){
		  double[][][] re = new double[sqVx.length][sqVx[0].length][sqVx[0][0].length];
		  for (int k = 0; k < sqVx.length  ; k++)
			  for (int j = 0 ; j < sqVx[0].length ; j++)
				  for ( int i = 0 ; i < sqVx[0][0].length  ; i++)
					  re[k][j][i] =(1.0 /Math.pow(deltaXY,2)) 
					  * Math.pow( 1.0 + sqVx[k][j][i] / Math.pow(deltaXY,2), -0.5 ) 
					  * Math.exp(  -  ( Math.sqrt(1.0 + sqVx[k][j][i] / Math.pow(deltaXY,2)) - 1.0)   );
		  return re;
	  }

	  /**
	   * This static method corresponds to 
	   * the option 7 of  of the w function,
	   * the penalty function 
	   * @param sqVx square of Vx input 3D data
	   * @param deltaXY double parameter
	   * @return The phi function of Vx square
	   */
	  public static double[][][] w_V_7(double[][][] sqVx, double deltaXY){
		  double[][][] re = new double[sqVx.length][sqVx[0].length][sqVx[0][0].length];
		  final double  pow = 0.25;
		  for (int k = 0; k < sqVx.length  ; k++)
			  for (int j = 0 ; j < sqVx[0].length ; j++)
				  for ( int i = 0 ; i < sqVx[0][0].length  ; i++)
					  re[k][j][i] =( (2.0*pow) / Math.pow(deltaXY,2)) * Math.pow(1.0 + sqVx[k][j][i] / Math.pow(deltaXY,2), pow-1.0 );
		  return re;
	  }


	  /**
	   * This static method dotMult realizes the multiplication 
	   * of the first 1D input data by the second 1D input data, 
	   * elements by elements
	   * @param a first 1D input data
	   * @param b the second 1D input data
	   * @return The multiplication product in a
	   */
	  public static double[] dotMult(double[] a, double[] b){
		  if (a.length != b.length )
			  IJ.error( "MMCal error Matrix shoult be same size." );
		  for (int k = 0; k < a.length  ; k++) 
			  a[k] *= b[k];
		  return a;
	  }

	  /**
	   * This static method dotMult realizes the multiplication 
	   * of the first 3D input data by the second 3D input data, 
	   * elements by elements
	   * @param a first 3D input data
	   * @param b the second 3D input data
	   * @return The multiplication product in a
	   */
	  public static double[][][] dotMult(double[][][] a, double[][][] b){
		  if (a.length != b.length || a[0].length != b[0].length || a[0][0].length != b[0][0].length)
			  IJ.error( "MMCal error Matrix shoult be same size." );
		  for (int k = 0; k < a.length  ; k++)
			  for (int j = 0 ; j < a[0].length ; j++)
				  for (int i = 0 ; i < a[0][0].length ; i++) 
					  a[k][j][i] *= b[k][j][i];
		  return a;
	  }

	  /**
	   * This static method dotMult realizes the multiplication 
	   * of the first 3D input data by the second 3D input data, 
	   * elements by elements
	   * @param a first 3D input data
	   * @param b the second 3D input data
	   * @return The multiplication product in a
	   */
	  public static double[][][] dotMult(double[][][] a, int[][][] b){
		  if (a.length != b.length || a[0].length != b[0].length || a[0][0].length != b[0][0].length)
			  IJ.error( "MMCal error Matrix shoult be same size." );
		  for (int k = 0; k < a.length  ; k++)
			  for (int j = 0 ; j < a[0].length ; j++)
				  for (int i = 0 ; i < a[0][0].length ; i++) 
					  a[k][j][i] *= b[k][j][i];
		  return a;
	  }

	  /**
	   * This static method  dotMultComplex realizes the multiplication 
	   * of the first 3D Complex input data by the second 3D Complex input data, 
	   * elements by elements 
	   * the input arrays must be of size slices by rows by 2*columns, 
	   * physical layout of the input data is as follows:
	   * <pre>
	   * a[k3][k2][2*k1] = Re[k3][k2][k1], a[k3][k2][2*k1+1] = Im[k3][k2][k1],
	   * 0&lt;=k3&lt;slices, 0&lt;=k2&lt;columns, 0&lt;=k1&lt;rows,
	   * </pre> 
	   * @param a first 3D complex input data
	   * @param b the second 3D complex input data
	   * @return The multiplication product in a
	   */
	  public static double[][][] dotMultComplex(double[][][] a, double[][][] b){
		  if (a.length != b.length || a[0].length != b[0].length || a[0][0].length != b[0][0].length)
			  IJ.error( "MMCal error Matrix shoult be same size." );
		  double[][][] ret = new double[a.length][a[0].length][a[0][0].length];
		  for (int k = 0; k < a.length  ; k++)
			  for (int j = 0 ; j < a[0].length ; j++)
				  for (int i = 0 ; i < a[0][0].length/2 ; i++) {
					  ret[k][j][2*i] = (a[k][j][2*i] *  b[k][j][2*i]) - (a[k][j][2*i+1] * b[k][j][2*i+1]);
					  ret[k][j][2*i+1] = (a[k][j][2*i] *  b[k][j][2*i+1]) + (a[k][j][2*i+1] * b[k][j][2*i]);
				  }
		  return ret;
	  }

	  /**
	   * This static method dotMult realizes the multiplication 
	   * of the first 3D input data by the value of <code>a</code>
	   * for each element
	   * @param x first 3D input data
	   * @param a the second 3D input data
	   * @return The multiplication product in x
	   */
	  public static double[][][] dotMult(double[][][] x, double a){
		  for (int k = 0; k < x.length  ; k++)
			  for (int j = 0 ; j < x[0].length ; j++)
				  for ( int i = 0 ; i < x[0][0].length  ; i++)
					  x[k][j][i] = a *  x[k][j][i];
		  return x;
	  }

	  /**
	   * This static method dotDiv realizes the division 
	   * of the first 3D input data by the second 3D input data, 
	   * elements by elements
	   * @param a first 3D input data
	   * @param b the second 3D input data
	   * @return The division result in a
	   */
	  public static double[][][] dotDiv(double[][][] a, double[][][] b){
		  if (a.length != b.length || a[0].length != b[0].length || a[0][0].length != b[0][0].length)
			  IJ.error( "MMCal error Matrix shoult be same size." );
		  for (int k = 0; k < a.length  ; k++)
			  for (int j = 0 ; j < a[0].length ; j++)
				  for (int i = 0 ; i < a[0][0].length ; i++) 
					  a[k][j][i] = a[k][j][i] / b[k][j][i];
		  return a;
	  }

	  /**
	   * This static method exp calculs the exponential
	   * of the  3D input data 
	   * @param x first 3D input data
	   * @return The exp result in x
	   */
	  public static double[][][] exp(double[][][] x){
		  for (int k = 0; k < x.length  ; k++)
			  for (int j = 0 ; j < x[0].length ; j++)
				  for ( int i = 0 ; i < x[0][0].length  ; i++)
					  x[k][j][i] = Math.exp(x[k][j][i]);
		  return x;
	  }

	  /**
	   * This static method add calculs the addition
	   * of the first 3D input data by the second 3D input data, 
	   * elements by elements
	   * @param a first 3D input data
	   * @param b the second 3D input data
	   * @return The addition result in a
	   */  
	  public static double[][][] add(double[][][] a, double[][][] b){
		  if (a.length != b.length || a[0].length != b[0].length || a[0][0].length != b[0][0].length)
			  IJ.error( "MMCal error Matrix shoult be same size." );
		  for (int k = 0; k < a.length  ; k++)
			  for (int j = 0 ; j < a[0].length ; j++)
				  for ( int i = 0 ; i < a[0][0].length  ; i++)
					  a[k][j][i] += b[k][j][i];
		  return a;
	  }

	  /**
	   * This static method minus calculs the subtraction
	   * of the first 3D input data by the second 3D input data, 
	   * elements by elements
	   * @param a first 3D input data
	   * @param b the second 3D input data
	   * @return The subtraction result in a
	   */
	  public static double[][][] minus(double[][][] a, double[][][] b){
		  if (a.length != b.length || a[0].length != b[0].length || a[0][0].length != b[0][0].length)
			  IJ.error( "MMCal error Matrix shoult be same size." );
		  for (int k = 0; k < a.length  ; k++)
			  for (int j = 0 ; j < a[0].length ; j++)
				  for ( int i = 0 ; i < a[0][0].length  ; i++)
					  a[k][j][i] -= b[k][j][i];
		  return a;
	  }   

	  /**
	   * This static method add calculs the addition
	   * of the  3D input data by the value  of <code>b</code>
	   * @param a first 3D input data
	   * @param b the double value to add
	   * @return The addition result in a
	   */  
	  public static double[][][] add(double[][][] a, double b){
		  for (int k = 0; k < a.length  ; k++)
			  for (int j = 0 ; j < a[0].length ; j++)
				  for ( int i = 0 ; i < a[0][0].length  ; i++)
					  a[k][j][i] += b;
		  return a;
	  }  

	  /**
	   * This static method add calculs the addition
	   * of the first 1D input data by the second 1D input data, 
	   * elements by elements
	   * @param a first 1D input data
	   * @param b the second 1D input data
	   * @return The addition result in a
	   */ 
	  public static double[] add(double[] a, double[] b){
		  for (int k = 0; k < a.length  ; k++)
			  a[k] += b[k];
		  return a;
	  } 

	  /**
	   * This static method minus calculs the subtraction
	   * of the first 1D input data by the second 1D input data, 
	   * elements by elements
	   * @param a first 1D input data
	   * @param b the second 1D input data
	   * @return The subtraction result in a
	   */
	  public static double[] minus(double[] a, double[] b){
		  for (int k = 0; k < a.length  ; k++)
			  a[k] -= b[k];
		  return a;
	  } 


	  /**
	   * This static method Voperator provides operator for Gradient in each direction
	   * The Outpout is composed of three arrays for each direction
	   * @param x 3D input data
	   * @return Array list contains Vvx, Vhx, Vtx Gradiant in each direction
	   */
	  public static ArrayList Voperator(double[][][] x){
		  ArrayList li = new ArrayList();
		  int Nx = x[0][0].length;
		  int Ny = x[0].length;
		  int Nz = x.length;
		  double[][][] Vvx = MMCal.copy(x);
		  double[][][] Vhx = MMCal.copy(x);
		  double[][][] Vtx = MMCal.copy(x);
		  for (int k = 0; k < Nz  ; k++)
			  for (int j = 0 ; j < Ny ; j++)
				  for ( int i = 0 ; i < Nx  ; i++)
				  {   
					  if(i != 0)
						  Vhx[k][j][i] = x[k][j][i] - x[k][j][i-1];
					  if(j != 0)
						  Vvx[k][j][i] = x[k][j][i] - x[k][j-1][i];
					  if(k != 0)
						  Vtx[k][j][i] = x[k][j][i] - x[k-1][j][i];
				  }    	
		  li.add(Vhx);
		  li.add(Vvx);
		  li.add(Vtx);
		  return li;
	  }


	  /**
	   * This static method Vvoperatoradj  provides operator for adjacent Gradient columns direction
	   * @param x 3D input data
	   * @return 3D output data, adjacent gradiant
	   */
	  public static double[][][] Vvoperatoradj(double[][][] x){
		  int Ny = x[0].length;
		  double[][][] Vvtx = MMCal.copy(x);    	
		  for (int k = 0; k < x.length  ; k++)
			  for (int j = 0 ; j < Ny ; j++)
				  for ( int i = 0 ; i < x[0][0].length  ; i++)
				  {   
					  if(j != Ny-1)
						  Vvtx[k][j][i] = x[k][j][i] - x[k][j+1][i];
				  }
		  return Vvtx;
	  }

	  /**
	   * This static method Vhoperatoradj  provides operator for adjacent Gradient rows direction
	   * @param x 3D input data
	   * @return 3D output data, adjacent gradiant
	   */
	  public static double[][][] Vhoperatoradj(double[][][] x){
		  int Nx = x[0][0].length;
		  double[][][] Vvtx = MMCal.copy(x);    	
		  for (int k = 0; k < x.length  ; k++)
			  for (int j = 0 ; j < x[0].length ; j++)
				  for ( int i = 0 ; i < Nx  ; i++)
				  {   
					  if(i != Nx-1)
						  Vvtx[k][j][i] = x[k][j][i] - x[k][j][i+1];
				  }
		  return Vvtx;
	  }

	  /**
	   * This static method Vtoperatoradj  provides operator for adjacent Gradient slices direction
	   * @param x 3D input data
	   * @return 3D output data, adjacent gradiant
	   */
	  public static double[][][] Vtoperatoradj(double[][][] x){
		  int Nz = x.length;
		  double[][][] Vttx = MMCal.copy(x);    	
		  for (int k = 0; k < Nz  ; k++)
			  for (int j = 0 ; j < x[0].length ; j++)
				  for ( int i = 0 ; i < x[0][0].length   ; i++)
				  {   
					  if(k != Nz-1)
						  Vttx[k][j][i] = x[k][j][i] - x[k+1][j][i];
				  }
		  return Vttx;
	  }


	  /**
	   * This static method norm calculs the norm 2 of the input data
	   * @param x 3D input data
	   * @return double value of the input norm two
	   */
	  public static double norm(double[] x){
		  double sum = 0;
		  for (int k = 0; k < x.length  ; k++)
			  sum = sum + Math.pow(Math.abs(x[k]) , 2);
		  sum = Math.sqrt(sum);

		  return sum;
	  } 


	  /**
	   * This static method norm calcul the norm "fro" otherwise norm 2  
	   * @param x 1D input data
	   * @param st name in String of the method's name "fro" implemented otherwise norm 2  
	   * @return double value of the input norm "fro" otherwise norm 2  
	   */
	  public static double norm(double[] x, String st){
		  double sum = 0;
		  if ( st.equalsIgnoreCase("fro") ){
			  for (int k = 0; k < x.length  ; k++)
				  sum = sum + Math.pow(Math.abs(x[k]) , 2);
			  sum = Math.sqrt(sum);
		  }
		  else {
			  for (int k = 0; k < x.length  ; k++)
				  sum = sum + Math.pow(Math.abs(x[k]) , 2);
			  sum = Math.sqrt(sum);

		  }
		  return sum;
	  }

	  /**
	   * This static method norm calcul the norm "fro" otherwise norm 2  
	   * @param x 3D input data
	   * @param st name in String of the method's name "fro" implemented otherwise norm 2  
	   * @return double value of the input norm "fro" otherwise norm 2  
	   */
	  public static double norm(double[][][] x, String st){
		  double sum = 0;
		  if ( st.equalsIgnoreCase("fro") ){
			  for (int k = 0; k < x.length  ; k++)
				  for (int j = 0 ; j < x[0].length ; j++)
					  for ( int i = 0 ; i < x[0][0].length   ; i++)
						  sum = sum + Math.pow(Math.abs(x[k][j][i]) , 2);
			  sum = Math.sqrt(sum);
		  }
		  else {
			  for (int k = 0; k < x.length  ; k++)
				  for (int j = 0 ; j < x[0].length ; j++)
					  for ( int i = 0 ; i < x[0][0].length   ; i++)
						  sum = sum + Math.pow(Math.abs(x[k][j][i]) , 2);
			  sum = Math.sqrt(sum);

		  }
		  return sum;
	  }

	  /**
	   * This Class implements methods to reduce dimension of the input data
	   * The reduction can concern different ways and is equivalent to the
	   * matlab squeeze function
	   * Squeeze sq = new MMCal.Squeeze();
	   * sq.cal(m);
	   * double  inter =  (double) sq.getSqueeze().get(0);
	   * or
	   * double[]  inter =  (double[]) sq.getSqueeze().get(0);
	   * or
	   * double[][]  inter =  (double[][]) sq.getSqueeze().get(0);
	   * or
	   * double[][][]  inter =  (double[][][]) sq.getSqueeze().get(0);
	   * Copyright:    Copyright (c) 2016
	   * Company:      Labex Archimède
	   *               AMU Aix Marseille Université
	   * @author       Dominique Benielli
	   * @author       Caroline Chaux-Moulin
	   * @author       Sandrine Anthoine
	   * @version 0.0 
	   */
	  public static class Squeeze {  	
		  private  double[][][] tab3d;
		  private  double[][] tab2d;
		  private  double[] tab1d;
		  private  double tab = Double.NaN;

		  public Squeeze(){
			  ;
		  }


		  /**
		   * This method gets one of the fourth possible results of Squeeze
		   * The First one is in 3D dimension
		   * The Second in in 2D dimension 
		   * The Third  in 1D dimension
		   * The Fourth is a number 
		   * Only one is in the output list
		   * @return the Firts element of list contains the data in reduced dimension
		   */
		  public  ArrayList getSqueeze(){
			  ArrayList li = new ArrayList();
			  if (this.tab3d != null){
				  li.add (this.tab3d);}
			  if (this.tab2d != null){
				  li.add (this.tab2d);}
			  if (this.tab1d != null)
				  li.add (this.tab1d);
			  if ( !Double.isNaN(this.tab)){
				  li.add (this.tab);}
			  return li;
		  }

		  /**
		   * This static method realizes the the squeeze of the 3D input data, 
		   * if some dimension are equal to one, the dimension to squeeze are specified
		   *  a == 0 && b == 0 dimensions to squezze rows and columns
		   *  a == 1 && b == 0 dimensions to squezze slices and columns
		   *  a == 0 && b == 1 dimensions to squezze slices and rows
		   * @param A 3D input data
		   * @param a 0 or 1 
		   * @param b 0 or 1 
		   * @return the squeezed output in 1D dimension
		   */
		  private static double[] rep(double[][][] A, int a, int b){
			  double[] ret = new double[0];
			  if ( a == 0 && b == 0){
				  ret = new double[A.length];
				  for(int k=0 ; k< A.length ; k++)
					  ret[k] = A[k][0][0];
			  }
			  if ( a == 1 && b == 0){
				  ret = new double[A[0][0].length];
				  for(int i=0 ; i< A[0][0].length ; i++)
					  ret[i] = A[0][0][i];
			  }
			  if ( a == 0 && b == 1){
				  ret = new double[A[0].length];
				  for(int j=0 ; j< A[0].length ; j++)
					  ret[j] = A[0][j][0];
			  }
			  return ret;
		  }

		  /**
		   * This static method realizes the the squeeze of the 3D input data, 
		   * if some dimension are equal to one, the dimension to squeeze are specified
		   *  a == 0 && b == 0 dimensions to squezze rows and columns
		   *  a == 1 && b == 0 dimensions to squezze slices and columns
		   *  a == 0 && b == 1 dimensions to squezze slices and rows
		   * @param A 3D input data
		   * @param a 0 or 1 
		   * @param b 0 or 1 
		   * @return the squeezed output in 1D dimension
		   */
		  private static double[][] rep(double[][][] A, int a){
			  double[][] ret = new double[0][0];     
			  switch (a){
			  case 0:
				  ret = new double[A[0].length][A.length] ;
				  break;
			  case 1:
				  ret = new double[A[0][0].length][A.length] ;
				  break;
			  case 2:
				  ret = new double[A[0].length][A[0][0].length] ;
				  break;        	
			  }
			  for(int k=0 ; k< A.length ; k++)
				  for (int j=0 ; j< A[0].length ; j++)
					  for (int i=0 ; i< A[0][0].length ; i++)
					  {
						  switch (a){
						  case 0:
							  ret[j][k] = A[k][j][0];
							  break;
						  case 1:
							  ret[i][k] = A[k][0][i];
							  break;
						  case 2:
							  ret[j][i] = A[0][j][i];
							  break;						   
						  }
					  }
			  return ret;
		  }

		  /**
		   * This static method realizes the the squeeze of the 2D input data, 
		   * if some dimension are equal to one, the dimension to squeeze are specified
		   *  a == 0  dimensions to squezze rows 
		   *  a == 1  dimensions to squezze columns
		   * @param A 2D input data
		   * @param a 0 or 1 
		   * @return the squeezed output in 1D dimension
		   */
		  private static double[] rep(double[][] A, int a){
			  double[] ret = new double[0];     
			  switch (a){
			  case 1:
				  ret = new double[A.length] ;
				  break;
			  case 0:
				  ret = new double[A[0].length] ;
				  break;      	
			  }
			  for (int j=0 ; j< A.length ; j++)
				  for (int i=0 ; i< A[0].length ; i++)
				  {
					  switch (a){
					  case 1:
						  ret[j] = A[j][0];
						  break;
					  case 0:
						  ret[i] = A[0][i];
						  break;
					  }
				  }
			  return ret;
		  }

		  /**
		   * This static method realizes the the squeeze of the 3D input data,
		   * use the 'getSqueeze()' method to to get the result
		   * @param A 3D input data
		   */
		  public void cal(double[][][] A){
			  // case 1 all != 1
			  if (A[0][0].length != 1 && A[0].length != 1 && A.length != 1){
				  this.tab3d = A;
				  this.tab2d =null;
				  this.tab1d =null;
				  this.tab = Double.NaN;
			  }

			  // case i == 1, j ang k != 1
			  if ((A[0][0].length == 1) && (A[0].length != 1) && (A.length != 1)){       		
				  this.tab2d = MMCal.Squeeze.rep( A, 0 );
				  this.tab3d =null;
				  this.tab1d =null;
				  this.tab = Double.NaN;   		    
			  }       	        	
			  // case j == 1, i ang k != 1
			  if ((A[0].length == 1) && (A[0][0].length != 1) && (A.length != 1)){     		
				  this.tab2d = MMCal.Squeeze.rep( A, 1 );
				  this.tab1d =null;
				  this.tab3d =null;
				  this.tab = Double.NaN;   		    
			  }
			  // case k == 1, i ang j != 1
			  if ((A[0][0].length != 1) && (A[0].length != 1) && (A.length == 1)){     		
				  this.tab2d = MMCal.Squeeze.rep( A, 2 );
				  this.tab1d =null;
				  this.tab3d =null;
				  this.tab = Double.NaN;   		    
			  }
			  // case i == 1 , j == 1 and k!= 1
			  if ((A[0][0].length == 1) && (A[0].length == 1) && (A.length != 1)){ 
				  this.tab1d = rep( A, 0, 0 );
				  this.tab2d =null;
				  this.tab3d =null;
				  this.tab = Double.NaN; 
			  }

			  // case k == 1 , j == 1 and i!= 1
			  if ((A[0][0].length != 1) && (A[0].length == 1) && (A.length == 1)){ 
				  this.tab1d = rep( A, 1,0 );
				  this.tab2d =null;
				  this.tab3d =null;
				  this.tab = Double.NaN; 
			  } 

			  // case k == 1 , i == 1 and j != 1
			  if ((A[0][0].length == 1) && (A[0].length != 1) && (A.length == 1)){ 
				  this.tab1d = rep( A, 0, 1 );
				  this.tab2d =null;
				  this.tab3d =null;
				  this.tab = Double.NaN; 
			  }
			  // case k == 1 , i == 1 and j 1= 1
			  if ((A[0][0].length == 1) && (A[0].length == 1) && (A.length == 1)){ 
				  this.tab1d = null;
				  this.tab2d =null;
				  this.tab3d =null;
				  this.tab = A[0][0][0]; 
			  }
		  }

		  /**
		   * This static method realizes the the squeeze of the 2D input data,
		   * use the 'getSqueeze()' method to to get the result
		   * @param A 2D input data
		   */
		  public void cal(double[][] A){

			  // case i == 1, j  != 1
			  if ((A[0].length == 1) && (A.length != 1)){       		
				  this.tab1d = rep( A, 1  );
				  this.tab3d =null;
				  this.tab2d =null;
				  this.tab = Double.NaN;   		    
			  }       	        	

			  //  j == 1 and i!= 1
			  if ( (A[0].length != 1) && (A.length == 1)){ 
				  this.tab1d = rep( A, 0 );
				  this.tab2d =null;
				  this.tab3d =null;
				  this.tab = Double.NaN; 
			  } 

			  // case j == 1 , i == 1 
			  if ((A[0].length == 1) && (A.length == 1)){ 
				  this.tab1d = null;
				  this.tab2d =null;
				  this.tab3d =null;
				  this.tab = A[0][0]; 
			  }
		  }  

		  /**
		   * This static method realizes the the squeeze of the 1D input data,
		   * use the 'getSqueeze()' method to to get the result
		   * @param A 1D input data
		   */
		  public void cal(double[] A){

			  // case i == 1
			  if ((A.length == 1)){       		
				  this.tab2d = null;
				  this.tab3d =null;
				  this.tab1d =null;
				  this.tab = A[0];   		    
			  } 
			  else{
				  this.tab1d = A;
				  this.tab2d =null;
				  this.tab3d =null;
				  this.tab = Double.NaN;
			  }

		  }  
	  }
}

