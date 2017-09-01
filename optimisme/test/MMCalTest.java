package optimisme.test;

import java.util.ArrayList;
import java.util.Arrays;
import optimisme.MMCal;
import optimisme.MMCal.Squeeze;

import org.junit.After;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;


public class MMCalTest {
	private static final double[][][] mat333 = {{{1.0,2.0,3.0},{1.0,2.0,3.0},{1.0,2.0,3.0}},
        {{4.0,5.0,6.0},{4.0,5.0,6.0},{4.0,5.0,6.0}},
        {{5.0,5.0,5.0},{5.0,5.0,5.0},{5.0,5.0,5.0}}};
	private static final double[][] mat32 = {{1.0,2.0,3.0},{1.0,2.0,3.0},{1.0,2.0,3.0},
        {4.0,4.0,4.0},{4.0,4.0,4.0},{4.0,4.0,4.0}};
	private static final double[][][] mat112 = {{{1.0}},{{1.0}}};
	private static final double[] mat13 = {1, 2, 3, 4};
	private static final double[][][] H = {{{0.01, 0.020, 0.01, 0.020, 0.03, 0.01},
		                           		{0.01, 0.020, 0.1, 0.02, 0.03, 0.01},
		                           		{0.01, 0.020, 0.01, 0.2, 0.03, 0.01}},
		                           {{0.01, 0.020, 0.01, 0.020, 0.03, 0.01},
		                           		{0.01, 1.0, 2.0, 3.0, 1.0, 0.020},
		                           		{0.01, 0.01, 0.2, 0.03, 0.01, 0.01}},
		                           {{0.01, 0.01, 0.020, 0.03, 0.01, 0.01},
		                           			{ 0.1, 1.0, 2.0, 4.0, 1.0, 0.020},
		                           			{0.1, 0.01, 0.2 ,0.03, 0.01, 0.01}},
		                           	{{0.01, 0.01, 0.020, 0.03, 0.01, 0.01},
		                           			{0.01, 0.1, 0.02, 0.03 ,0.01, 0.020},
		                           			{0.1, 0.01, 0.2, 0.03, 0.01, 0.01}}};

	@Before
	public void setUpBefore() throws Exception {
		
	}

	@After
	public void tearDownAfter() throws Exception {
	}
	
	@Test
	public void testConcatAll333() {
		double[] m = MMCal.concatAll(mat333);
		double[] mexpected = {1.0,1.0 ,1.0 , 2.0, 2.0,2.0, 3.0, 3.0 ,3.0,
		        4.0,4.0, 4.0, 5.0 , 5.0 ,5.0 ,6.0,6.0,6.0,
		        5.0,5.0,5.0, 5.0,5.0,5.0,5.0,5.0,5.0};
		Assert.assertArrayEquals(mexpected, m,10E-2);
	}
	
	@Test
	public void testLTE333() {
		int[][][] m = MMCal.lte(3.0, mat333);
		int[][][] mexpected = {{{1,1,1},{1,1,1},{1,1,1}},
		        {{0,0,0},{0,0,0},{0,0,0}},
		        {{0,0,0},{0,0,0},{0,0,0}}};;
		Assert.assertArrayEquals(mexpected, m);
	}
	
	@Test
	public void testGT333() {
		int[][][] m = MMCal.gt(3.0, mat333);
		int[][][] mexpected = {{{0,0,0},{0,0,0},{0,0,0}},
		        {{1,1,1},{1,1,1},{1,1,1}},
		        {{1,1,1},{1,1,1},{1,1,1}}};;
		Assert.assertArrayEquals(mexpected, m);
	}
	
	@Test
	public void testGTE333() {
		int[][][] m = MMCal.gte(3.0, mat333);
		int[][][] mexpected = {{{0,0,1},{0,0,1},{0,0,1}},
		        {{1,1,1},{1,1,1},{1,1,1}},
		        {{1,1,1},{1,1,1},{1,1,1}}};;
		Assert.assertArrayEquals(mexpected, m);
	}

	@Test
	public void testVvoperatoradj333() {
		double[][][] Vvtx = MMCal.Vvoperatoradj(H);
		double[][][] Vvtxexpected = {{{ 
			0.00000 ,  0.00000  ,-0.09000  , 0.00000 ,  0.00000 ,  0.00000},
			{   0.00000 ,  0.00000  , 0.09000 , -0.18000 ,  0.00000  , 0.00000},
			{  0.01000  , 0.02000  , 0.01000  , 0.20000 ,  0.03000 ,  0.01000}},

			{{   0.00000 , -0.98000  ,-1.99000  ,-2.98000  ,-0.97000 , -0.01000},
				{  0.00000,   0.99000 ,  1.80000  , 2.97000  , 0.99000  , 0.01000},
				{  0.01000  , 0.01000 ,  0.20000 ,  0.03000 ,  0.01000  , 0.01000}},

				{{   -0.09000 , -0.99000 , -1.98000,  -3.97000 , -0.99000  ,-0.01000},
					{   0.00000 ,  0.99000  , 1.80000  , 3.97000  , 0.99000 ,  0.01000},
					{  0.10000 ,  0.01000,   0.20000  , 0.03000  , 0.01000  , 0.01000}},

					{ { 0.00000 , -0.09000 ,  0.00000 ,  0.00000  , 0.00000 , -0.01000},
						{  -0.09000 ,  0.09000 , -0.18000 ,  0.00000 ,  0.00000 ,  0.01000},
						{ 0.10000  , 0.01000 ,  0.20000  , 0.03000 ,  0.01000  , 0.01000}}};

		for (int k =0 ; k<  Vvtx.length;k++){
			for (int j =0 ; j<  Vvtx[0].length; j++){
				Assert.assertArrayEquals( Vvtxexpected[k][j],  Vvtx[k][j], 10-6);
			}
		}
	}
	
	@Test
	public void testconverttoFloat2d(){
		float[][] factual = MMCal.converttoFloat2d(mat32);
		float[][] fexpected = {{(float)1.0,(float)2.0,(float)3.0},
				{(float)1.0,(float)2.0,(float)3.0},{(float)1.0,(float)2.0,(float)3.0},
		        {(float)4.0,(float) 4.0,(float)4.0},{(float)4.0,(float)4.0,(float)4.0},
		        {(float)4.0,(float)4.0,(float)4.0}};
		Assert.assertArrayEquals(fexpected,factual);
	}
	
	@Test
	public void testconverttoFloat1d(){
		float[] factual = MMCal.converttoFloat1d(mat32);
		float[] fexpected = {(float)1.0,(float)2.0,(float)3.0,
				(float)1.0,(float)2.0,(float)3.0,(float)1.0,(float)2.0,(float)3.0,
		        (float)4.0,(float) 4.0,(float)4.0,(float)4.0,(float)4.0,(float)4.0,
		        (float)4.0,(float)4.0,(float)4.0};
		Assert.assertArrayEquals(fexpected,factual,(float)10e-3);
	}
	
	@Test
	public void testgetWithFilter(){
		
		int[][][] filter = {{{1,1,1},{0,0,0},{1,1,0}},
		        {{0,0,0},{0,0,0},{0,1,1}},
		        {{0,0,1},{1,1,1},{0,0,0}}};
		double[] factual = MMCal.getWithFilter(mat333, filter);
		double[] fexpected = {1.0,2.0,3.0,1.0,2.0, 5.0,6.0,5.0,5.0,5.0,5.0};
		Assert.assertArrayEquals(fexpected,factual,10e-2);
		int[][][] filter2 = {{{0,0,0},{0,0,0},{0,0,0}},
		        {{0,0,0},{0,0,0},{0,0,0}},
		        {{0,0,0},{0,0,0},{0,0,0}}};
		factual = MMCal.getWithFilter(mat333, filter2);
		double[] fexpected2 = {};
		Assert.assertArrayEquals(fexpected2,factual,10e-2);
	}
	
	@Test
	public void testsetValue2d(){
		
		double[][] input = MMCal.copy(mat32);
		input = MMCal.setValue(input, 6.0);
		
		double[][] fexpected ={{6.0,6.0,6.0},{6.0,6.0,6.0},{6.0,6.0,6.0},
		        {6.0,6.0,6.0},{6.0,6.0,6.0},{6.0,6.0,6.0}};

		Assert.assertArrayEquals(fexpected,input);


	}
	
	@Test
	public void testsetValue1d(){
		
		double[] input = Arrays.copyOf(mat13,mat13.length);
		input = MMCal.setValue(input, 6.0);
		
		double[] fexpected  = {6.0, 6.0, 6, 6};

		Assert.assertArrayEquals(fexpected,input,10e-2);


	}
	
	@Test
	public void testVhoperatoradj333() {
		double[][][] Vvtx = MMCal.Vhoperatoradj(H);
		double[][][] Vvtxexpected = {{{ 
			-0.010000  , 0.010000 , -0.010000 ,-0.010000  , 0.020000  , 0.010000},
			{ -0.010000 , -0.080000  , 0.080000 , -0.010000  , 0.020000  , 0.010000},
			{ -0.010000 ,  0.010000 , -0.190000 ,  0.170000 ,  0.020000  , 0.010000}},

			{{  -0.01000,   0.01000 , -0.01000  ,-0.01000  , 0.02000 ,  0.01000},
			{  -0.99000 , -1.00000  ,-1.00000  , 2.00000  , 0.98000 ,  0.02000},
			{   0.00000 , -0.19000  , 0.17000 ,  0.02000  , 0.00000 ,  0.01000}},

			{{   0.00000 , -0.01000  ,-0.01000  , 0.02000  , 0.00000 ,  0.01000},
			{   -0.90000 , -1.00000 , -2.00000  , 3.00000 ,  0.98000 ,  0.02000},
			{    0.09000 , -0.19000 ,  0.17000  , 0.02000 ,  0.00000,   0.01000}},

			{{  0.00000 , -0.01000 , -0.01000 ,  0.02000 ,  0.00000 ,  0.01000},
			{  -0.09000 ,  0.08000 , -0.01000 ,  0.02000 , -0.01000 ,  0.02000},
			{  0.09000 , -0.19000 ,  0.17000  , 0.02000 ,  0.00000  , 0.01000}}};

		for (int k =0 ; k<  Vvtx.length;k++){
			for (int j =0 ; j<  Vvtx[0].length; j++){
			
				Assert.assertArrayEquals( Vvtxexpected[k][j],  Vvtx[k][j], 10-6);
			}
		}
	}

	@Test
	public void testVtoperatoradj333() {
		double[][][] Vttx = MMCal.Vtoperatoradj(H);
		double[][][] Vttxexpected = {{{ 
			   0.00000 ,  0.00000 ,  0.00000 ,  0.00000 ,  0.00000 ,  0.00000},
			   {  0.00000 , -0.98000  ,-1.90000 , -2.98000 , -0.97000 , -0.01000},
			   {  0.00000  , 0.01000 , -0.19000 ,  0.17000  , 0.02000 ,  0.00000	}},

			{{    0.00000 ,  0.01000 , -0.01000  ,-0.01000  , 0.02000  , 0.00000},
			{  -0.09000 ,  0.00000 ,  0.00000 , -1.00000 ,  0.00000  , 0.00000},
			 {  -0.09000 ,  0.00000 ,  0.00000 ,  0.00000 ,  0.00000 ,  0.00000}},

			{{   0.00000,   0.00000 ,  0.00000  , 0.00000  , 0.00000 ,  0.00000},
			 {0.09000  , 0.90000 ,  1.98000  , 3.97000  , 0.99000  , 0.00000},
			{0.00000  , 0.00000 , 0.00000 ,  0.00000  , 0.00000 ,  0.00000  }},

			{{    0.010000 ,  0.010000  , 0.020000 , 0.030000  , 0.010000  , 0.010000},
			{ 0.010000  , 0.100000 ,  0.020000  , 0.030000  , 0.010000 ,  0.020000},
			{ 0.100000 ,  0.010000  , 0.200000 ,  0.030000 ,  0.010000  , 0.010000}}};

		for (int k =0 ; k<  Vttx.length;k++){
			for (int j =0 ; j<  Vttx[0].length; j++){
			
				Assert.assertArrayEquals( Vttxexpected[k][j],  Vttx[k][j], 10-6);
			}
		}
	}
	
	
	@Test
	public void testPhi_sqV_1() {
		double deltaXY =  0.01;
		double[][][] sqVx = MMCal.square( MMCal.copy(H)); 
		double[][][] phi_sqV_1expected = {{  
			{ 0.39347 ,  0.86466 ,  0.39347 ,  0.86466 ,  0.98889  , 0.39347},
				{ 0.39347 , 0.86466 ,  1.00000 ,  0.86466 ,  0.98889  , 0.39347},
				{  0.39347 ,  0.86466 ,  0.39347 ,  1.00000 ,  0.98889 ,  0.39347}},
			{ 
				{  0.39347 ,  0.86466 ,  0.39347 ,  0.86466  , 0.98889,   0.39347},
				{  0.39347 ,  1.00000 ,  1.00000 ,  1.00000  , 1.00000 ,  0.86466},
				{  0.39347 ,  0.39347  , 1.00000 ,  0.98889  , 0.39347 ,  0.39347}},
			{
				{ 	0.39347 ,  0.39347 ,  0.86466 ,  0.98889 ,  0.39347 ,  0.39347},
				{   1.00000 ,  1.00000 ,  1.00000 ,  1.00000 ,  1.00000 ,  0.86466},
				{   1.00000 ,  0.39347 ,  1.00000 ,  0.98889 ,  0.39347  , 0.39347}},
			{ 
				{   0.39347 ,  0.39347 ,  0.86466  , 0.98889  , 0.39347 ,  0.39347},
				{   0.39347 ,  1.00000 ,  0.86466 ,  0.98889  , 0.39347  , 0.86466},
				{  1.00000  , 0.39347  , 1.00000 ,  0.98889  , 0.39347  , 0.39347}}};
			
						
		double[][][] result = MMCal.phi_sqV_1(sqVx, deltaXY);
		for (int k =0 ; k< result.length;k++)
			for (int j =0 ; j< result[0].length; j++)
					Assert.assertArrayEquals(phi_sqV_1expected[k][j],  result[k][j],10-2);		
		}
	

	@Test
	public void testPhi_sqV_2() {
		double deltaXY =  0.01;
		double[][][] sqVx = MMCal.square( MMCal.copy(H)); 
		double[][][] phi_sqV_2expected = {{  
			{ 0.33333  , 0.66667  , 0.33333 ,  0.66667 ,  0.81818 ,  0.33333},
			{ 0.33333 ,  0.66667 ,  0.98039 ,  0.66667 ,  0.81818,   0.33333},
			{ 0.33333 ,  0.66667 ,  0.33333 , 0.99502 ,0.81818  , 0.33333}},
		{
			{ 0.33333  , 0.66667 ,  0.33333  , 0.66667  , 0.81818 ,  0.33333},
			{ 0.33333  , 0.99980 ,  0.99995  , 0.99998 ,  0.99980 ,  0.66667},
			{ 0.33333  , 0.33333 ,  0.99502  , 0.81818 , 0.33333 ,  0.33333}},
		{
			{ 0.33333 ,  0.33333 ,  0.66667 , 0.81818  , 0.33333 ,  0.33333},
			{ 0.98039 ,  0.99980 ,  0.99995 ,  0.99999 ,  0.99980,  0.66667},
			{ 0.98039  , 0.33333 ,  0.99502 ,  0.81818 , 0.33333 , 0.33333}},
	  {
			{ 0.33333  , 0.33333 ,  0.66667 ,  0.81818 ,  0.33333  , 0.33333},
			{ 0.33333 ,  0.98039 ,  0.66667  , 0.81818  ,0.33333 ,  0.66667},
			{ 0.98039 ,  0.33333 ,  0.99502 ,  0.81818 , 0.33333 ,  0.33333}}};
					   					
		double[][][] result = MMCal.phi_sqV_2(sqVx, deltaXY);
		for (int k =0 ; k< result.length;k++)
			for (int j =0 ; j< result[0].length; j++)
					Assert.assertArrayEquals(phi_sqV_2expected[k][j],  result[k][j],10-2);		
		}
	
	@Test
	public void testPhi_sqV_3() {
		double deltaXY =  0.01;
		double[][][] sqVx = MMCal.square( MMCal.copy(H)); 
		double[][][] phi_sqV_3expected = {{  
			       {0.69315 ,  1.60944 ,  0.69315  , 1.60944   ,2.30259,   0.69315},
				   {0.69315 ,  1.60944  , 4.61512,   1.60944 ,  2.30259,   0.69315},
				   {0.69315 ,  1.60944  , 0.69315 ,  5.99396 ,  2.30259  , 0.69315}},
				{
				  { 0.69315  ,  1.60944  ,  0.69315   , 1.60944  ,  2.30259  ,  0.69315},
				  { 0.69315  ,  9.21044 ,  10.59666  , 11.40758 ,  9.21044  ,  1.60944},
				  { 0.69315 ,   0.69315 ,   5.99396  ,  2.30259  ,  0.69315 ,   0.69315}},
				{
				  {0.69315   , 0.69315  ,  1.60944 ,   2.30259   , 0.69315  , 0.69315},
				  { 4.61512 ,   9.21044 ,  10.59666  , 11.98294  ,  9.21044 ,  1.60944},
				  { 4.61512 ,   0.69315 ,   5.99396  ,  2.30259  ,  0.69315  ,  0.69315}},
				{
				  { 0.69315 ,  0.69315 ,  1.60944 ,  2.30259 ,  0.69315 ,  0.69315},
				  { 0.69315 ,  4.61512  , 1.60944 ,  2.30259 ,0.69315 , 1.60944},
				  { 4.61512  , 0.69315 ,  5.99396 ,  2.30259  , 0.69315  ,0.69315}}};
				  
		double[][][] result = MMCal.phi_sqV_3(sqVx, deltaXY);
		for (int k =0 ; k< result.length;k++)
			for (int j =0 ; j< result[0].length; j++)
					Assert.assertArrayEquals(phi_sqV_3expected[k][j],  result[k][j],10-2);		
		}
	@Test
	public void testPhi_sqV_4() {
		double deltaXY =  0.01;
		double[][][] sqVx = MMCal.square( MMCal.copy(H)); 
		double[][][] phi_sqV_4expected = {{  
				{   0.41421   , 1.23607   , 0.41421 ,  1.23607   , 2.16228  ,  0.41421},
				{    0.41421  ,  1.23607  ,  9.04988 ,   1.23607  ,  2.16228  ,  0.41421},
				{    0.41421   , 1.23607  ,  0.41421,  19.02498   , 2.16228  ,  0.41421}},
			{
				{     0.41421  ,   1.23607  ,   0.41421  ,   1.23607  ,  2.16228 ,    0.41421},
				{   0.41421    ,99.00500 ,  199.00250  , 299.00167 ,   99.00500  ,   1.23607},
				{   0.41421   ,  0.41421 ,   19.02498  ,   2.16228 ,    0.41421   ,  0.41421}},
			{
				{  0.41421   ,  0.41421  ,  1.23607  ,   2.16228  ,   0.41421  ,   0.41421},
				{  9.04988   , 99.00500 ,  199.00250 ,  399.00125 ,   99.00500 ,    1.23607},
				{  9.04988   ,  0.41421  ,  19.02498  ,   2.16228   ,  0.41421  ,   0.41421}},
			{
				{ 0.41421  ,  0.41421   , 1.23607 ,   2.16228  ,  0.41421  ,  0.41421},
				{   0.41421 ,   9.04988  ,  1.23607  ,  2.16228 ,   0.41421  , 1.23607},
				{ 9.04988   , 0.41421 ,  19.02498   , 2.16228 ,   0.41421 ,   0.41421}}};
							
		double[][][] result = MMCal.phi_sqV_4(sqVx, deltaXY);
		for (int k =0 ; k< result.length;k++)
			for (int j =0 ; j< result[0].length; j++)
					Assert.assertArrayEquals(phi_sqV_4expected[k][j],  result[k][j],10-2);		
		}
	
	@Test
	public void testPhi_sqV_5() {
		double deltaXY =  0.01;
		double[][][] sqVx = MMCal.square( MMCal.copy(H)); 
		double[][][] phi_sqV_5expected = MMCal.dotMult(MMCal.copy(sqVx), 0.5);
							
		double[][][] result = MMCal.phi_sqV_5(sqVx, deltaXY);
		for (int k =0 ; k< result.length;k++)
			for (int j =0 ; j< result[0].length; j++)
					Assert.assertArrayEquals(phi_sqV_5expected[k][j],  result[k][j],10-2);		
		}
	
	@Test
	public void testPhi_sqV_6() {
		double deltaXY =  0.01;
		double[][][] sqVx = MMCal.square( MMCal.copy(H)); 
		double[][][] phi_sqV_6expected = {{  
			{ 0.33914  , 0.70948 ,  0.33914 ,  0.70948  , 0.88494 ,  0.33914},
			{ 0.33914 ,  0.70948 ,  0.99988 ,  0.70948  , 0.88494 ,  0.33914},
			{ 0.33914  , 0.70948 ,  0.33914  , 1.00000 ,  0.88494 ,  0.33914}},
			{
			{ 0.33914 ,  0.70948 ,  0.33914  , 0.70948  , 0.88494 ,  0.33914},
			{ 0.33914 ,  1.00000 ,  1.00000 ,  1.00000 ,  1.00000 ,  0.70948},
			{ 0.33914  , 0.33914 ,  1.00000  , 0.88494  , 0.33914 ,  0.33914}},
			{
			{ 0.33914 ,  0.33914 ,  0.70948  , 0.88494 ,  0.33914 ,  0.33914},
			{ 0.99988  , 1.00000 ,  1.00000  , 1.00000 ,  1.00000 ,  0.70948},
			{ 0.99988 ,  0.33914  , 1.00000  , 0.88494 ,  0.33914 ,  0.33914}},
			{
			{ 0.33914 ,  0.33914  , 0.70948,   0.88494 ,  0.33914 ,  0.33914},
			{ 0.33914 ,  0.99988 , 0.70948  , 0.88494  , 0.33914 ,  0.70948},
			{ 0.99988 ,  0.33914  , 1.00000  , 0.88494,   0.33914 , 0.33914}}};
			
			
								
		double[][][] result = MMCal.phi_sqV_6(sqVx, deltaXY);
		for (int k =0 ; k< result.length;k++)
			for (int j =0 ; j< result[0].length; j++)
					Assert.assertArrayEquals(phi_sqV_6expected[k][j],  result[k][j],10-2);		
		}
	
	@Test
	public void testPhi_sqV_7() {
		double deltaXY =  0.01;
		double[][][] sqVx = MMCal.square( MMCal.copy(H)); 
		double[][][] phi_sqV_7expected = {{  			
			{  0.18921  , 0.49535 ,  0.18921 ,  0.49535 ,  0.77828 ,  0.18921},
			{  0.18921  , 0.49535 ,  2.17015 ,  0.49535 ,  0.77828 , 0.18921},
			{  0.18921  , 0.49535 ,  0.18921 ,  3.47493,   0.77828  , 0.18921}},
			{   
			{  0.18921 ,   0.49535  ,  0.18921 ,   0.49535  ,  0.77828  ,  0.18921},
			{  0.18921  ,  9.00025 ,  13.14222 ,  16.32056  ,  9.00025  ,  0.49535},
			{  0.18921 ,   0.18921 ,   3.47493 ,  0.77828   , 0.18921  , 0.18921}},
			{   
			{  0.18921  ,  0.18921  ,  0.49535  ,  0.77828  , 0.18921 ,  0.18921},
			{  2.17015 ,   9.00025 , 13.14222 , 19.00003   , 9.00025  ,  0.49535},
			{  2.17015  ,  0.18921 ,   3.47493 ,   0.77828  ,  0.18921   , 0.18921}},
			{   
			{   0.18921,   0.18921  , 0.49535 , 0.77828 ,  0.18921,   0.18921},
			{  0.18921 ,  2.17015 ,  0.49535 , 0.77828 ,  0.18921 ,  0.49535},
			{  2.17015 ,  0.18921  , 3.47493 ,  0.77828 ,  0.18921 ,  0.18921}}};
			   
							
		double[][][] result = MMCal.phi_sqV_7(sqVx, deltaXY);
		for (int k =0 ; k< result.length;k++)
			for (int j =0 ; j< result[0].length; j++)
					Assert.assertArrayEquals(phi_sqV_7expected[k][j],  result[k][j],10-2);		
		}
	
	@Test
	public void testw_V_1() {
		double deltaXY =  0.01;
		double[][][] sqVx = MMCal.square( MMCal.copy(H)); 
		double[][][] w_V_1expected = {{  
		{ 6.0653e+03 ,  1.3534e+03  , 6.0653e+03  , 1.3534e+03  , 1.1109e+02 ,  6.0653e+03},
		{   6.0653e+03 ,  1.3534e+03 ,  1.9287e-18 ,  1.3534e+03 ,  1.1109e+02 ,  6.0653e+03},
		{ 6.0653e+03,   1.3534e+03  , 6.0653e+03 ,  1.3839e-83  , 1.1109e+02  , 6.0653e+03}},
	{
		{ 6.0653e+03 ,  1.3534e+03 ,  6.0653e+03 ,  1.3534e+03 ,  1.1109e+02  , 6.0653e+03},
		{ 6.0653e+03 ,  0.0000e+00 ,  0.0000e+00  , 0.0000e+00 ,  0.0000e+00 ,  1.3534e+03},
		{ 6.0653e+03,   6.0653e+03,   1.3839e-83 ,  1.1109e+02 ,  6.0653e+03 ,  6.0653e+03}},
	{
		{ 6.0653e+03 ,  6.0653e+03 ,  1.3534e+03 ,  1.1109e+02,   6.0653e+03  , 6.0653e+03},
		{   1.9287e-18 ,  0.0000e+00  , 0.0000e+00 ,  0.0000e+00 ,  0.0000e+00  , 1.3534e+03},
		{ 1.9287e-18  , 6.0653e+03 ,  1.3839e-83  , 1.1109e+02  , 6.0653e+03  , 6.0653e+03}},
	{
		{ 6.0653e+03 ,  6.0653e+03 ,  1.3534e+03  , 1.1109e+02 ,  6.0653e+03,   6.0653e+03},
		{ 6.0653e+03 ,  1.9287e-18 ,  1.3534e+03,  1.1109e+02,   6.0653e+03,  1.3534e+03},
		{ 1.9287e-18 ,  6.0653e+03 ,  1.3839e-83, 1.1109e+02 ,  6.0653e+03 ,  6.0653e+03}}};
			
		double[][][] result = MMCal.w_V_1(sqVx, deltaXY);
		for (int k =0 ; k< result.length;k++)
			for (int j =0 ; j< result[0].length; j++)
					Assert.assertArrayEquals(w_V_1expected[k][j],  result[k][j],10-2);		
		}
	
	@Test
	public void testw_V_2() {
		double deltaXY =  0.01;
		double[][][] sqVx = MMCal.square( MMCal.copy(H)); 
		double[][][] w_V_2expected = {{
				{	4.4444e+03  , 1.1111e+03  , 4.4444e+03 ,  1.1111e+03 ,  3.3058e+02 ,  4.4444e+03},
				{   4.4444e+03  , 1.1111e+03 ,  3.8447e+00 ,  1.1111e+03  , 3.3058e+02 ,  4.4444e+03},
				{   4.4444e+03  , 1.1111e+03 ,  4.4444e+03 ,  2.4752e-01 ,  3.3058e+02 ,  4.4444e+03}},
			{
				{   4.4444e+03  , 1.1111e+03 ,  4.4444e+03  , 1.1111e+03 ,  3.3058e+02 ,  4.4444e+03},
				{	4.4444e+03  , 3.9984e-04 ,  2.4998e-05 ,  4.9381e-06 ,  3.9984e-04  , 1.1111e+03},
				{	4.4444e+03  , 4.4444e+03 ,  2.4752e-01 ,  3.3058e+02 ,  4.4444e+03 ,  4.4444e+03}},

			{
				{	4.4444e+03  , 4.4444e+03  , 1.1111e+03  , 3.3058e+02 ,  4.4444e+03 ,  4.4444e+03},
				{	3.8447e+00  , 3.9984e-04 ,  2.4998e-05  , 1.5625e-06 ,  3.9984e-04 ,  1.1111e+03},
				{   3.8447e+00  , 4.4444e+03  , 2.4752e-01 ,  3.3058e+02 ,  4.4444e+03  , 4.4444e+03}},
			{
				{   4.4444e+03  , 4.4444e+03  , 1.1111e+03 ,  3.3058e+02  , 4.4444e+03 ,  4.4444e+03},
				{   4.4444e+03  , 3.8447e+00 ,  1.1111e+03 ,  3.3058e+02 ,  4.4444e+03 ,  1.1111e+03},
				{   3.8447e+00  , 4.4444e+03,   2.4752e-01 ,  3.3058e+02 ,  4.4444e+03  , 4.4444e+03}}};

			
		double[][][] result = MMCal.w_V_2(sqVx, deltaXY);
		for (int k =0 ; k< result.length;k++)
			for (int j =0 ; j< result[0].length; j++)
					Assert.assertArrayEquals(w_V_2expected[k][j],  result[k][j],10-2);		
		}
	
	@Test
	public void testw_V_3() {
		double deltaXY =  0.01;
		double[][][] sqVx = MMCal.square( MMCal.copy(H)); 
		double[][][] w_V_3expected = {{		
			{   1.0000e+04 ,  4.0000e+03 ,  1.0000e+04 ,  4.0000e+03 ,  2.0000e+03 ,  1.0000e+04},
			{  1.0000e+04 ,  4.0000e+03 ,  1.9802e+02 ,  4.0000e+03  , 2.0000e+03  , 1.0000e+04},
			{  1.0000e+04 ,  4.0000e+03 ,  1.0000e+04 ,  4.9875e+01 ,  2.0000e+03 ,  1.0000e+04}},
			{
			{ 1.0000e+04 ,  4.0000e+03  , 1.0000e+04 ,  4.0000e+03 ,  2.0000e+03  , 1.0000e+04},
			{ 1.0000e+04 ,  1.9998e+00 ,  4.9999e-01 ,  2.2222e-01 ,  1.9998e+00  , 4.0000e+03},
			{ 1.0000e+04 ,  1.0000e+04  , 4.9875e+01 ,  2.0000e+03 ,  1.0000e+04  , 1.0000e+04}},
			{
			{ 1.0000e+04,   1.0000e+04  , 4.0000e+03 ,  2.0000e+03 ,  1.0000e+04  , 1.0000e+04},
			{ 1.9802e+02 ,  1.9998e+00 ,  4.9999e-01 ,  1.2500e-01 ,  1.9998e+00 ,  4.0000e+03},
			{ 1.9802e+02  , 1.0000e+04 , 4.9875e+01  , 2.0000e+03  , 1.0000e+04  , 1.0000e+04}},
			{ 
			{   1.0000e+04 ,  1.0000e+04 ,  4.0000e+03 ,  2.0000e+03  , 1.0000e+04 , 1.0000e+04},
			{   1.0000e+04 ,  1.9802e+02  , 4.0000e+03 ,  2.0000e+03 ,  1.0000e+04 ,  4.0000e+03},
			{   1.9802e+02  , 1.0000e+04 ,  4.9875e+01 ,  2.0000e+03 , 1.0000e+04 ,  1.0000e+04}}};

			
		double[][][] result = MMCal.w_V_3(sqVx, deltaXY);
		for (int k =0 ; k< result.length;k++)
			for (int j =0 ; j< result[0].length; j++)
					Assert.assertArrayEquals(w_V_3expected[k][j],  result[k][j],10-2);		
		}
	
	@Test
	public void testw_V_4() {
		double deltaXY =  0.01;
		double[][][] sqVx = MMCal.square( MMCal.copy(H)); 
		double[][][] w_V_4expected = {{
			{7071.07 ,  4472.14  , 7071.07 ,  4472.14 ,  3162.28  , 7071.07},
			{7071.07 ,  4472.14  ,  995.04  , 4472.14  , 3162.28 ,  7071.07},
			{7071.07  , 4472.14  , 7071.07  ,  499.38 ,  3162.28 ,  7071.07}},
			{
			{7071.068  , 4472.136  , 7071.068  , 4472.136 ,  3162.278  , 7071.068},
			{7071.068  ,   99.995  ,   49.999  ,   33.333  ,   99.995  , 4472.136},
			{ 7071.068 ,  7071.068  ,  499.376 ,  3162.278 ,  7071.068 ,  7071.068}},
			{
			{7071.068  , 7071.068 , 4472.136 ,  3162.278 ,  7071.068 ,  7071.068},
			{995.037  ,   99.995  ,   49.999  ,   25.000  ,   99.995 ,  4472.136},
			{995.037  , 7071.068  ,  499.376 ,  3162.278 ,  7071.068 ,  7071.068}},
			{
			{7071.07  , 7071.07 ,  4472.14  , 3162.28  , 7071.07 ,  7071.07},
			{7071.07 ,   995.04 ,  4472.14 ,  3162.28  , 7071.07 ,  4472.14},
			{995.04  , 7071.07 ,   499.38  , 3162.28 ,  7071.07  , 7071.07}}};
			
		double[][][] result = MMCal.w_V_4(sqVx, deltaXY);
		for (int k =0 ; k< result.length;k++)
			for (int j =0 ; j< result[0].length; j++)
					Assert.assertArrayEquals(w_V_4expected[k][j],  result[k][j],10-2);		
		}
	
	
	@Test
	public void testw_V_5() {
		double deltaXY =  0.01;
		double[][][] sqVx = MMCal.square( MMCal.copy(H)); 
		double[][][] w_V_5expected = new double[H.length][H[0].length][H[0][0].length];
		MMCal.add(w_V_5expected , 1.0);   
		double[][][] result = MMCal.w_V_5(sqVx, deltaXY);
		for (int k =0 ; k< result.length;k++)
			for (int j =0 ; j< result[0].length; j++)
					Assert.assertArrayEquals(w_V_5expected[k][j],  result[k][j],10-2);
			
		}
	
	@Test
	public void testw_V_6() {
		double deltaXY =  0.01;
		double[][][] sqVx = MMCal.square( MMCal.copy(H)); 
		double[][][] w_V_6expected = {{
	   {4.6730e+03 ,  1.2993e+03  , 4.6730e+03 ,  1.2993e+03  , 3.6386e+02  , 4.6730e+03},
	   {4.6730e+03 ,  1.2993e+03 ,  1.1682e-01 ,  1.2993e+03 ,  3.6386e+02  , 4.6730e+03},
	   {4.6730e+03 ,  1.2993e+03  , 4.6730e+03  , 2.7289e-06 ,  3.6386e+02  , 4.6730e+03}},
		 {
	    {4.6730e+03  ,  1.2993e+03  ,  4.6730e+03  ,  1.2993e+03  ,  3.6386e+02  ,  4.6730e+03},
	    {4.6730e+03  ,  1.0061e-41  ,  1.8762e-85  , 4.6570e-129   , 1.0061e-41  ,  1.2993e+03},
	    {4.6730e+03  ,  4.6730e+03 ,   2.7289e-06  ,  3.6386e+02   , 4.6730e+03  ,  4.6730e+03}},
	    {
	    {4.6730e+03  ,  4.6730e+03 ,   1.2993e+03  ,  3.6386e+02 ,   4.6730e+03 ,   4.6730e+03},
	    {1.1682e-01  ,  1.0061e-41  ,  1.8762e-85 ,  1.2999e-172 ,   1.0061e-41 ,   1.2993e+03},
	    {1.1682e-01  ,  4.6730e+03 ,   2.7289e-06  ,  3.6386e+02 ,   4.6730e+03 ,   4.6730e+03}},
	    {
	   {4.6730e+03 ,  4.6730e+03 ,  1.2993e+03 ,  3.6386e+02  , 4.6730e+03 ,  4.6730e+03},
	   {4.6730e+03 ,  1.1682e-01 ,  1.2993e+03  , 3.6386e+02  , 4.6730e+03 ,  1.2993e+03},
	   {1.1682e-01 ,  4.6730e+03 ,  2.7289e-06 ,  3.6386e+02  , 4.6730e+03,   4.6730e+03}}};
	   
	   
	   
		double[][][] result = MMCal.w_V_6(sqVx, deltaXY);
		for (int k =0 ; k< result.length;k++)
			for (int j =0 ; j< result[0].length; j++)
					Assert.assertArrayEquals(w_V_6expected[k][j],  result[k][j],10-2);
			
		}
	
	@Test
	public void testw_V_7() {
		double deltaXY =  0.01;
		double[][][] sqVx = MMCal.square( MMCal.copy(H)); 
		double[][][] w_V_7expected = 
		{{
		{2973.018 , 1495.349  , 2973.018  , 1495.349  ,  889.140 ,  2973.018},  
	    {2973.018  , 1495.349 ,   156.938  , 1495.349  ,  889.140 ,  2973.018},
		{2973.018  , 1495.349 ,  2973.018  ,   55.797 ,   889.140  , 2973.018}},
		{
		 {2.9730e+03  , 1.4953e+03  , 2.9730e+03 ,  1.4953e+03 ,  8.8914e+02  , 2.9730e+03},
		 {2.9730e+03  , 4.9996e+00 ,  1.7677e+00 ,  9.6224e-01 ,  4.9996e+00  , 1.4953e+03},
		 {2.9730e+03 , 2.9730e+03 ,  5.5797e+01 ,  8.8914e+02  , 2.9730e+03   ,2.9730e+03} },
		{
		   {2.9730e+03 ,  2.9730e+03  , 1.4953e+03  , 8.8914e+02  , 2.9730e+03 ,  2.9730e+03},
		   { 1.5694e+02 ,  4.9996e+00 ,  1.7677e+00 ,  6.2500e-01 ,  4.9996e+00  , 1.4953e+03},
		   {1.5694e+02 ,  2.9730e+03 , 5.5797e+01 ,  8.8914e+02 ,  2.9730e+03  , 2.9730e+03}},
		{
		   {2973.018 ,  2973.018  , 1495.349  ,  889.140 ,  2973.018 ,  2973.018},
		   {2973.018 ,   156.938  , 1495.349   , 889.140 ,  2973.018 ,  1495.349},
		   {156.938  , 2973.018   ,  55.797  ,  889.140  , 2973.018 ,  2973.018}}};
	   
		double[][][] result = MMCal.w_V_7(sqVx, deltaXY);
		for (int k =0 ; k< result.length;k++)
			for (int j =0 ; j< result[0].length; j++)
					Assert.assertArrayEquals(w_V_7expected[k][j],  result[k][j],10-2);
			
		}
	
	@Test
	public void testVoperator333() {


		ArrayList li = MMCal.Voperator(H);
		double[][][] Vvx = (double[][][]) li.get(0);
		double[][][]  Vhx = (double[][][]) li.get(1);
		double[][][]  Vtx = (double[][][]) li.get(2);
		double[][][] Vvxexpected = {{{  0.01000 ,  0.02000,   0.01000 ,  0.02000  , 0.03000 ,  0.01000},
			{ 0.00000 ,  0.00000  , 0.09000  , 0.00000  , 0.00000  , 0.00000},
			{  0.00000  , 0.00000  ,-0.09000 ,  0.18000   ,0.00000 ,  0.00000}},

			{   {0.01000 ,  0.02000 ,  0.01000 ,  0.02000 ,  0.03000  , 0.01000},
				{   0.00000 ,  0.98000 ,  1.99000  , 2.98000 ,  0.97000 ,  0.01000},
				{	   0.00000 , -0.99000 , -1.80000 , -2.97000 , -0.99000 , -0.01000}},

				{ { 0.01000 ,  0.01000 ,  0.02000 ,  0.03000 ,  0.01000 ,  0.01000},
					{   0.09000 ,  0.99000  , 1.98000  , 3.97000 ,  0.99000  , 0.01000},
					{		   0.00000 , -0.99000 , -1.80000 , -3.97000 , -0.99000  ,-0.01000}},

					{  {  0.01000 ,  0.01000  , 0.02000 ,  0.03000 ,  0.01000  , 0.01000},
						{   0.00000,   0.09000 ,  0.00000 ,  0.00000 ,  0.00000  , 0.01000},
						{ 0.09000 , -0.09000,   0.18000,   0.00000,   0.00000,  -0.01000} }};

		double[][][] Vhxexpected = {
				{{   0.010000  , 0.010000 , -0.010000 ,  0.010000,   0.010000 , -0.020000},
					{   0.010000  , 0.010000 ,  0.080000  ,-0.080000 ,  0.010000 , -0.020000},
					{ 0.010000 ,  0.010000 , -0.010000  , 0.190000 , -0.170000 , -0.020000}},

					{	{   0.01000 ,  0.01000  ,-0.01000 ,  0.01000 ,  0.01000 , -0.02000},
						{   0.01000 ,  0.99000  , 1.00000 ,  1.00000 , -2.00000 , -0.98000},
						{ 0.01000  , 0.00000,   0.19000 , -0.17000 , -0.02000 ,  0.00000}},

						{{   0.01000 ,  0.00000  , 0.01000  , 0.01000 , -0.02000 ,  0.00000},
							{   0.10000 ,  0.90000  , 1.00000 ,  2.00000 , -3.00000 , -0.98000},
							{  0.10000 , -0.09000 ,  0.19000 , -0.17000 , -0.02000 ,  0.00000}},

							{{  0.01000 ,  0.00000  , 0.01000 ,  0.01000  ,-0.02000 ,  0.00000},
								{   0.01000 ,  0.09000 , -0.08000 , 0.01000 , -0.02000 ,  0.01000},
								{ 0.10000 , -0.09000  , 0.19000 , -0.17000 , -0.02000 ,  0.00000}}};

		double[][][] Vtxexpected = {
				{{   0.010000  , 0.020000 ,  0.010000  , 0.020000 ,  0.030000 ,  0.010000},
					{   0.010000 , 0.020000  , 0.100000  , 0.020000,   0.030000,   0.010000},
					{ 0.010000  , 0.020000 ,  0.010000 ,  0.200000 ,  0.030000 ,  0.010000}},

					{{   0.00000  , 0.00000  , 0.00000 ,  0.00000 ,  0.00000,   0.00000},
						{	0.00000 ,  0.98000  , 1.90000 ,  2.98000  , 0.97000  , 0.01000},
						{	0.00000,  -0.01000 ,  0.19000 , -0.17000  ,-0.02000  , 0.00000}},

						{{   0.00000 , -0.01000 ,  0.01000 ,  0.01000 , -0.02000  , 0.00000},
							{   0.09000 ,  0.00000  , 0.00000 ,  1.00000  , 0.00000 ,  0.00000},
							{   0.09000 ,  0.00000 ,  0.00000 ,  0.00000 ,  0.00000 ,  0.00000}},

							{{  0.00000   ,0.00000  , 0.00000  , 0.00000  , 0.00000  , 0.00000},
								{   -0.09000,  -0.90000 , -1.98000  ,-3.97000 , -0.99000  , 0.00000},
								{  0.00000 ,  0.00000 ,  0.00000  , 0.00000 ,  0.00000 ,  0.00000}}};


		//	System.out.print(String.valueOf(Vvxexpected.length)+ " "+String.valueOf(Vvx.length));
		for (int k =0 ; k< Vvx.length;k++){
			for (int j =0 ; j< Vvx[0].length; j++){
				Assert.assertArrayEquals(Vvxexpected[k][j], Vvx[k][j], 10-6);
				Assert.assertArrayEquals(Vhxexpected[k][j], Vhx[k][j], 10-6);
				Assert.assertArrayEquals(Vtxexpected[k][j], Vtx[k][j], 10-6);
			}
		}
		/*	double[][][] mexpected = {{{1.0, 2.0, 3.0}},
                                   {{4.0, 5.0, 6.0}},
                                   {{5.0, 5.0, 5.0}}};*/

		//		Assert.assertArrayEquals(mexpected, m);
		//fail("Not yet implemented");
	} 

	
	@Test
	public void testExp333() {
		double[][][] expected =  
		        {{{  2.7183 ,   7.3891,   20.0855},
		        	{ 2.7183  ,  7.3891 ,  20.0855},
		        	{  2.7183  ,  7.3891  , 20.0855}},
		{{ 54.598  , 148.413 ,  403.429},
		{    54.598  , 148.413  , 403.429},
			{     54.598 ,  148.413 ,  403.429}},	
		{{   148.41 ,  148.41 ,  148.41},
		{   148.41  , 148.41 ,  148.41},
		{   148.41  , 148.41  , 148.41}}};
		double[][][] m = MMCal.exp(MMCal.copy(mat333));
    	for (int k=0 ; k< m.length ; k++)
    		for(int j=0 ; j< m[0].length ; j++)
    			Assert.assertArrayEquals(expected[k][j], m[k][j],10e-3);
	}
	
	@Test
	public void testNormElse333() {
		double expected =  22.3159136044;
		double m = MMCal.norm(mat333,"fro");
		Assert.assertEquals(expected, m, 10e-6);
		m = MMCal.norm(mat333,"pow");
		Assert.assertEquals(expected, m, 10e-6);
	}
	
	@Test
	public void testNorm13() {
		double expected =  5.477225575;
		double m = MMCal.norm(mat13);
		Assert.assertEquals(expected, m, 10e-6);
	}
	
	@Test
	public void testNormFro13() {
		double expected =  5.477225575;
		double m = MMCal.norm(mat13,"fro");
		Assert.assertEquals(expected, m, 10e-6);
	}
	
	@Test
	public void testNormElse13() {
		double expected =  5.477225575;
		double m = MMCal.norm(mat13,"pow");
		Assert.assertEquals(expected, m, 10e-6);
	}
	
	@Test
	public void testSquare13() {
		double[] m = MMCal.square(Arrays.copyOf(mat13,mat13.length));
		double[] mexpected  =  {1, 4, 9, 16};
		Assert.assertArrayEquals(mexpected, m, 10e-2);
	}
	
	@Test
	public void testSquare333() {
		double[][][] m = MMCal.square(MMCal.copy(mat333));
		double[][][] mexpected  = {{{1.0,4.0,9.0},{1.0,4.0,9.0},{1.0,4.0,9.0}},
		        {{16.0,25.0,36.0},{16.0,25.0,36.0},{16.0,25.0,36.0}},
		        {{25.0,25.0,25.0},{25.0,25.0,25.0},{25.0,25.0,25.0}}};
		Assert.assertArrayEquals(mexpected, m);
	}
	
	@Test
	public void testSqrt333() {	
		double[][][] mexpected  = {{{1.0,4.0,9.0},{1.0,4.0,9.0},{1.0,4.0,9.0}},
		        {{16.0,25.0,36.0},{16.0,25.0,36.0},{16.0,25.0,36.0}},
		        {{25.0,25.0,25.0},{25.0,25.0,25.0},{25.0,25.0,25.0}}};
		double[][][] m = MMCal.sqrt(mexpected);
		Assert.assertArrayEquals(m, mat333);
	}
	
	@Test
	public void testPow333() {
		double[][][] m = MMCal.pow(MMCal.copy(mat333),3);
		double[][][] mexpected  = 
			{{{1.0,8.0,27.0},{1.0,8.0,27.0},{1.0,8.0,27.0}},
		        {{64.0,125.0,216.0},{64.0,125.0,216.0},{64.0,125.0,216.0}},
		        {{125.0,125.0,125.0},{125.0,125.0,125.0},{125.0,125.0,125.0}}};
		Assert.assertArrayEquals(mexpected, m);
	}
	
	@Test
	public void testLog333() {
		double[][][] m = MMCal.log(MMCal.copy(mat333));
		double[][][] mexpected  = 
			{{{ 0.00000 ,  0.69315  , 1.09861},
			{0.00000 ,  0.69315 ,  1.09861},
			{ 0.00000 ,  0.69315  , 1.09861}},
			{{1.3863 ,  1.6094 ,  1.7918},
			{1.3863  , 1.6094  , 1.7918},
			{1.3863 , 1.6094 ,  1.7918}},  
			{{ 1.6094,   1.6094,   1.6094},
			{1.6094  , 1.6094 ,  1.6094},
			{1.6094  , 1.6094  , 1.6094}}};
    	for (int k=0 ; k< m.length ; k++)
    		for(int j=0 ; j< m[0].length ; j++)
    				Assert.assertArrayEquals(mexpected[k][j], m[k][j], 10e-5);
	}
	
	@Test
	public void testdotMultDouble13() {
		double[] m = MMCal.dotMult(Arrays.copyOf(mat13,mat13.length), mat13);
		double[] mexpected  =  {1, 4, 9, 16};
		Assert.assertArrayEquals(mexpected, m, 10e-2);
	}
	
	@Test(expected=Exception.class)
	public void testdotMultFaillure13() { 
    double[] mm = new double[mat13.length-2];
	MMCal.dotMult(Arrays.copyOf(mat13,mat13.length), mm);
	}
	
	@Test
	public void testdotMult333() {
		double[][][] m = MMCal.dotMult(MMCal.copy(mat333), 2);
		double[][][] mexpected  = {{{2.0,4.0,6.0},{2.0,4.0,6.0},{2.0,4.0,6.0}},
			        {{8.0,10.0,12.0},{8.0,10.0,12.0},{8.0,10.0,12.0}},
			        {{10.0,10.0,10.0},{10.0,10.0,10.0},{10.0,10.0,10.0}}};
		Assert.assertArrayEquals(mexpected, m);
	}
	
	@Test
	public void testdotMultDouble333() {
		double[][][] m = MMCal.dotMult(MMCal.copy(mat333), mat333);
		double[][][] mexpected  = {{{1.0,4.0,9.0},{1.0,4.0,9.0},{1.0,4.0,9.0}},
		        {{16.0,25.0,36.0},{16.0,25.0,36.0},{16.0,25.0,36.0}},
		        {{25.0,25.0,25.0},{25.0,25.0,25.0},{25.0,25.0,25.0}}};
		Assert.assertArrayEquals(mexpected, m);
	}
	
	@Test
	public void testdotMultInt333() { 
    int[][][] mm = new int[mat333.length][mat333[0].length][mat333[0][0].length];
    for( int k =0; k < mm.length ; k++)
    	for( int j=0 ; j< mm[0].length ; j++)
    			Arrays.fill(mm[k][j], 1);
	double[][][] mexpected = MMCal.dotMult(MMCal.copy(mat333), mm);
	Assert.assertArrayEquals(mexpected, mat333);
	}
	
	@Test(expected=Exception.class)
	public void testdotMultIntFaillure333() { 
    int[][][] mm = new int[mat333.length][mat333[0].length-1][mat333[0][0].length];
	MMCal.dotMult(MMCal.copy(mat333), mm);
	}
	
	@Test(expected=Exception.class)
	public void testdotMultDoubleFaillure333() { 
    double[][][] mm = new double[mat333.length][mat333[0].length-1][mat333[0][0].length];
	MMCal.dotMult(MMCal.copy(mat333), mm);
	}
	
	@Test
	public void testdotDivDouble333() {
		double[][][] m = MMCal.dotDiv(MMCal.copy(mat333), mat333);
		double[][][] mexpected  = {{{1.0,1.0,1.0},{1.0,1.0,1.0},{1.0,1.0,1.0}},
		        {{1.0,1.0,1.0},{1.0,1.0,1.0},{1.0,1.0,1.0}},
		        {{1.0,1.0,1.0},{1.0,1.0,1.0},{1.0,1.0,1.0}}};
		Assert.assertArrayEquals(mexpected, m);
	}

	@Test(expected=Exception.class)
	public void testdotDivDoubleFaillure333() { 
    double[][][] mm = new double[mat333.length][mat333[0].length-1][mat333[0][0].length];
	MMCal.dotDiv(MMCal.copy(mat333), mm);
	}
			
	@Test
	public void testMax13() {
		double[] m = MMCal.max(mat13, 3);
		double[] mexpected = {3, 3, 3, 4};
		Assert.assertArrayEquals(mexpected, m, 10e-4);
	}
	
	@Test
	public void testMaxInt13() {
		int[] mat13 = {1, 2, 3, 4};
		int[] m = MMCal.max(mat13, 3);
		int[] mexpected = {3, 3, 3, 4};
		Assert.assertArrayEquals(mexpected, m);
	}
	
	@Test
	public void testMax32() {
		double[][] m = MMCal.max(mat32, 3);
		double[][] mexpected = {{3.0,3.0,3.0},{3.0,3.0,3.0},{3.0,3.0,3.0},
		        {4.0,4.0,4.0},{4.0,4.0,4.0},{4.0,4.0,4.0}};
		Assert.assertArrayEquals(mexpected, m);
	}
	
	@Test
	public void testMax32String() {
		double[][] m = MMCal.max(mat32,"[]",0);
		double[][] mexpected = {{3.0},{3.0},{3.0},
		        {4.0},{4.0},{4.0}};
		Assert.assertArrayEquals(mexpected, m);
		
		 m = MMCal.max(mat32,"[]",1);
		 double[][] mexpected2 = {{4.0,4.0,4.0}};
		Assert.assertArrayEquals(mexpected2, m);
	
	}
	
	@Test
	public void testMax33() {
		double[][][] m = MMCal.max(mat333, 3);
		double[][][] mexpected = {{{3.0,3.0,3.0},{3.0,3.0,3.0},{3.0,3.0,3.0}},
		        {{4.0,5.0,6.0},{4.0,5.0,6.0},{4.0,5.0,6.0}},
		        {{5.0,5.0,5.0},{5.0,5.0,5.0},{5.0,5.0,5.0}}};
		Assert.assertArrayEquals(mexpected, m);
	}
	
	@Test
	public void testMax333() {
		double[][][] m = MMCal.max(mat333);
		double[][][] mexpected = {{{1.0,2.0,3.0}},
                {{4.0,5.0,6.0}},
                {{5.0,5.0,5.0}}};
		Assert.assertArrayEquals(mexpected, m);
	}
	
	@Test
	public void testMax0() {
		double[][][] mat = {{{1,3,6,7},{4,2,9,1}}};
		double[][][] m = MMCal.max(mat);
		double[][][] mexpected = {{{4,3,9,7}}};
		Assert.assertArrayEquals(mexpected, m);
	}

	@Test
	public void testMax1() {
		double[][][] mat = {{{1}},{{6}}};
		double[][][] m = MMCal.max(mat);
		double[][][] mexpected = {{{6}}};
		Assert.assertArrayEquals(mexpected, m);
	}
	
	@Test
	public void testMax2() {
		double[][][] mat = {{{1,3,6,7}}};
		double[][][] m = MMCal.max(mat);
		double[][][] mexpected = {{{7}}};
		Assert.assertArrayEquals(mexpected, m);
	}
	
	@Test
	public void testMaxElse() {
		double[][][] mat = {{{7}}};
		double[][][] m = MMCal.max(mat);
		double[][][] mexpected = {{{7}}};
		Assert.assertArrayEquals(mexpected, m);
	}
	
	@Test
	public void testmaxWithIndex1(){
		double[][][] m = (double[][][]) MMCal.maxWithIndex(mat333,1).get(0);
		int[][][] id = (int[][][]) MMCal.maxWithIndex(mat333,1).get(1);
		double[][][] mexpected =
			 {{{5.0,5.0,6.0},{5.0,5.0,6.0},{5.0,5.0,6.0}}};
		int[][][] idexpected =  {{{2,1,1},{2,1,1},{2,1,1}}};

		Assert.assertArrayEquals( mexpected,m);
		Assert.assertArrayEquals( idexpected,id);
	}	
	
	@Test
	public void testmaxWithIndex2(){
		double[][][] m = (double[][][]) MMCal.maxWithIndex(mat333,2).get(0);
		int[][][] id = (int[][][]) MMCal.maxWithIndex(mat333,2).get(1);
		double[][][] mexpected ={{{3.0},{3.0},{3.0}},
		        {{6.0},{6.0},{6.0}},
		        {{5.0},{5.0},{5.0}}};

		int[][][] idexpected =  {{{2},{2},{2}},{{2},{2},{2}},{{0},{0},{0}}};

		Assert.assertArrayEquals( mexpected,m);
		Assert.assertArrayEquals( idexpected,id);
	}	
	
	@Test
	public void testMin13() {	
	double[]  m = MMCal.min(mat13,3);
	double[]  mexpected =  {1, 2, 3, 3};
	Assert.assertArrayEquals(mexpected, m, 10e-6);
	}
	
	@Test
	public void testMinInt13() {	
	int[] mm = { 2 , 3, 6,5,7,10};
	int[]  m = MMCal.min(mm,5);
	int[]  mexpected =  { 2 , 3, 5,5,5,5};
	Assert.assertArrayEquals(mexpected, m);
	}
	
	@Test
	public void testAdd333() {
		double[][][] m = MMCal.add(MMCal.copy(mat333), 2);
		double[][][] mexpected  = {{{3.0,4.0,5.0},{3.0,4.0,5.0},{3.0,4.0,5.0}},
	        {{6.0,7.0,8.0},{6.0,7.0,8.0},{6.0,7.0,8.0}},
	        {{7.0,7.0,7.0},{7.0,7.0,7.0},{7.0,7.0,7.0}}};
		Assert.assertArrayEquals(mexpected, m);
	}	
	
	@Test(expected=Exception.class)
	public void testAddDoubleFaillure333() { 
    double[][][] mm = new double[mat333.length][mat333[0].length-1][mat333[0][0].length];
	MMCal.add(MMCal.copy(mat333), mm);
	}
	
	@Test
	public void testAddDouble333() { 
		double[][][] m = MMCal.add(MMCal.copy(mat333), mat333);

		double[][][] mexpected = {{{2.0,4.0,6.0},{2.0,4.0,6.0},{2.0,4.0,6.0}},
	        {{8.0,10.0,12.0},{8.0,10.0,12.0},{8.0,10.0,12.0}},
	        {{10.0,10.0,10.0},{10.0,10.0,10.0},{10.0,10.0,10.0}}};
		Assert.assertArrayEquals(mexpected, m);
	}
	
	@Test
	public void testAdd13() {	
	double[]  m = MMCal.add(Arrays.copyOf(mat13, mat13.length),mat13);
	double[]  mexpected = {2, 4, 6, 8};
	Assert.assertArrayEquals(mexpected, m, 10e-6);
	}	
	
	@Test
	public void testMinus13() {	
	double[]  m = MMCal.minus(Arrays.copyOf(mat13, mat13.length),mat13);
	double[]  mexpected = {0, 0, 0, 0};
	Assert.assertArrayEquals(mexpected, m, 10e-6);
	}	
	
	@Test
	public void testMinusDouble333() { 
		double[][][] m = MMCal.minus(MMCal.copy(mat333), mat333);
		double[][][] mexpected = {{{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}},
	        {{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}},
	        {{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}}};
		Assert.assertArrayEquals(mexpected, m);
	}
	
	@Test(expected=Exception.class)
	public void testMinusDoubleFaillure333() { 
    double[][][] mm = new double[mat333.length][mat333[0].length-1][mat333[0][0].length];
	MMCal.minus(MMCal.copy(mat333), mm);
	}
	
	
	@Test
	public void testSum12() {
		double[][] input = {{1,3,4}};
		double[][] m = MMCal.sum(input);
		double[][] mexpected = {{8}};
		Assert.assertArrayEquals(mexpected, m);
		double[][] input2 = {{1},{3},{4}};
		m = MMCal.sum(input2);
		double[][] mexpected2 = {{8}};
		Assert.assertArrayEquals(mexpected2, m);
		double[][] input3 = {{3}};
		m = MMCal.sum(input3);
		double[][] mexpected3 = {{3}};
		Assert.assertArrayEquals(mexpected3, m);
		
		
	}	
	
	@Test
	public void test0Sum32() {
		double[][] m = MMCal.sum(mat32,1);
		double[][] mexpected = {{6.0},{6.0},{6.0},{12.0},{12.0},{12.0}};

		Assert.assertArrayEquals(mexpected, m);
	}	

	@Test
	public void test1Sum32() {
		double[][] m = MMCal.sum(mat32,0);
		double[][] mexpected = {{15.0, 18.0, 21.0}};

		Assert.assertArrayEquals(mexpected, m);
	}
	
	@Test
	public void testSum333() {
		double[][][] m = MMCal.sum(mat333);
		double[][][] mexpected = {{{3.0,6.0,9.0}},
                                   {{12.0,15.0,18.0}},
                                   {{15.0,15.0,15.0}}};

		Assert.assertArrayEquals(mexpected, m);
	}
	
	@Test
	public void test2Sum333() {
		double[][][] m = MMCal.sum(mat333,2);
		double[][][] mexpected = {{{10.0,12.0,14.0}
                ,{10.0,12.0,14.0}
                ,{10.0,12.0,14.0}}};	
		Assert.assertArrayEquals(mexpected, m);
	}
	
	@Test
	public void testSum112() {

		double[][][] m = MMCal.sum(mat112);
		double[][][] mexpected ={{{2.0}}};	
		Assert.assertArrayEquals(mexpected, m);
		double[][][] m2 = {{{1}},{{2}}};
		m = MMCal.sum(m2);
		double[][][] mexpected2 = {{{3}}};
		Assert.assertArrayEquals(mexpected2, m);
		double[][][] m3 = {{{3}}};
		m = MMCal.sum(m3);
		double[][][] mexpected3 = {{{3}}};
		Assert.assertArrayEquals(mexpected3, m);
	}
	
	@Test
	public void testSum32() {

		double[][][] m = MMCal.sum(mat112);
		double[][][] mexpected ={{{2.0}}};	
		Assert.assertArrayEquals(mexpected, m);
		double[][][] m2 = {{{1}},{{2}}};
		m = MMCal.sum(m2);
		double[][][] mexpected2 = {{{3}}};
		Assert.assertArrayEquals(mexpected2, m);
		double[][][] m3 = {{{3}}};
		m = MMCal.sum(m3);
		double[][][] mexpected3 = {{{3}}};
		Assert.assertArrayEquals(mexpected3, m);
	}	
	

	@Test
	public void test1Sum333() {
		double[][][] m = MMCal.sum(mat333,1);
		double[][][] mexpected = {{{6.0},{6.0},{6.0}}
                ,{{15.0},{15.0},{15.0}}
                ,{{15.0},{15.0},{15.0}}};	
		Assert.assertArrayEquals(mexpected, m);
	}

	@Test
	public void test0Sum333() {
		double[][][] m = MMCal.sum(mat333,0);
		double[][][] mexpected = {{{3.0, 6.0 ,9.0}}
                ,{{12.0 ,15.0 ,18.0}}
                ,{{15.0 ,15.0 ,15.0}}};	
		Assert.assertArrayEquals(mexpected, m);
	}
	
	@Test
	public void test2Squezze() {
		double[][][] m =  {{{1.0,2.0,3.0},{1.0,2.0,3.0},{1.0,2.0,3.0},
		        {4.0,4.0,4.0},{4.0,4.0,4.0},{4.0,4.0,4.0}}};
		Squeeze sq = new MMCal.Squeeze();
		sq.cal(m);
		double [][] inter =  (double[][]) sq.getSqueeze().get(0);
		Assert.assertArrayEquals(inter, mat32);
	}

	@Test
	public void test1Squezze() {
		double[][][] m =  {{{1.0,2.0,3.0},{1.0,2.0,3.0},{1.0,2.0,3.0},
		        {4.0,4.0,4.0},{4.0,4.0,4.0},{4.0,4.0,4.0}}};
		Squeeze sq = new MMCal.Squeeze();
		sq.cal(m);
		double [][] inter =  (double[][]) sq.getSqueeze().get(0);
		Assert.assertArrayEquals(inter, mat32);
	}

	@Test
	public void test33Squezze() {
		double[][][] m =  {{{1.0}},{{2.0}},{{3.0}}};
		double[] mexpected ={1.0,2.0,3.0};
		Squeeze sq = new MMCal.Squeeze();
		sq.cal(m);
		double [] inter =  (double[]) sq.getSqueeze().get(0);
		Assert.assertArrayEquals(inter,mexpected,10-2);
	}
	
	@Test
	public void test32Squezze() {
		double[][][] m =  {{{1.0},{2.0},{3.0}}};
		double[] mexpected ={1.0,2.0,3.0};
		Squeeze sq = new MMCal.Squeeze();
		sq.cal(m);
		double [] inter =  (double[]) sq.getSqueeze().get(0);
		Assert.assertArrayEquals(inter,mexpected,10-2);
	}
	
	@Test
	public void test22Squezze() {
		double[][] m =  {{1.0},{2.0},{3.0}};
		double[] mexpected ={1.0,2.0,3.0};
		Squeeze sq = new MMCal.Squeeze();
		sq.cal(m);
		double [] inter =  (double[]) sq.getSqueeze().get(0);
		Assert.assertArrayEquals(inter,mexpected,10-2);
	}
	
	@Test
	public void test31Squezze() {
		double[][][] m =  {{{1.0,2.0,3.0}}};
		double[] mexpected ={1.0,2.0,3.0};
		Squeeze sq = new MMCal.Squeeze();
		sq.cal(m);
		double [] inter =  (double[]) sq.getSqueeze().get(0);
		Assert.assertArrayEquals(inter,mexpected,10-2);
	}
	
	@Test
	public void test11Squezze() {
		double[][][] m =  {{{1.0}}};
		double mexpected =1.0;
		Squeeze sq = new MMCal.Squeeze();
		sq.cal(m);
		double  inter =  (double) sq.getSqueeze().get(0);
		Assert.assertEquals(inter,mexpected,10e-2);
	}
	
	@Test
	public void test21Squezze() {
		double[][] m =  {{1.0}};
		double mexpected =1.0;
		Squeeze sq = new MMCal.Squeeze();
		sq.cal(m);
		double  inter =  (double) sq.getSqueeze().get(0);
		Assert.assertEquals(inter,mexpected,10e-2);
	}
	
	@Test
	public void test210Squezze() {
		double[] m =  {1.0};
		double mexpected =1.0;
		Squeeze sq = new MMCal.Squeeze();
		sq.cal(m);
		double  inter =  (double) sq.getSqueeze().get(0);
		Assert.assertEquals(inter,mexpected,10e-2);
	}
	
	@Test
	public void testNo2Squezze() {
		double[] m =  {1.0,2.0};
		double[] mexpected ={1.0,2.0};
		Squeeze sq = new MMCal.Squeeze();
		sq.cal(m);
		double[]  inter =  (double[]) sq.getSqueeze().get(0);
		Assert.assertArrayEquals(inter,mexpected,10e-2);
	}
	
	@Test
	public void testNo3Squezze() {

		Squeeze sq = new MMCal.Squeeze();
		sq.cal(mat333);
		double [][][] inter =  (double[][][]) sq.getSqueeze().get(0);
		Assert.assertArrayEquals(inter, mat333);
	}
	
	@Test
	public void test0Squezze() {
		double[][][] m =  {{{1.0},{2.0},{3.0}},{{1.0},{2.0},{3.0}},{{1.0},{2.0},{3.0}},
		        {{4.0},{4.0},{4.0}},{{4.0},{4.0},{4.0}},{{4.0},{4.0},{4.0}}};
		Squeeze sq = new MMCal.Squeeze();
		sq.cal(m);
		double [][] inter =  (double[][]) sq.getSqueeze().get(0);
		Assert.assertArrayEquals(inter, MMCal.nonConjugateTranspose(mat32));
	}
	
	@Test
	public void testnonConjugateTranspose() {
		double[][] m =  {{1.0, 1.0, 1.0, 4.0, 4.0, 4.0},
				           {2.0, 2.0, 2.0, 4.0, 4.0, 4.0},
				           {3.0, 3.0, 3.0, 4.0, 4.0, 4.0}};

		Assert.assertArrayEquals(m, MMCal.nonConjugateTranspose(mat32));
	}
	
	@Test
	public void x0Test(){
    	double[][][] Hb = H.clone();   
    	double maxhi =  MMCal.max(MMCal.concatAll(H))[0];
    	Hb = MMCal.setValueWithFilter(0.0, H, MMCal.lt(0.5*maxhi,H)); 
		Squeeze sq = new MMCal.Squeeze();
        sq.cal(MMCal.sum(Hb,2));
        
        double[][] inter2= (double[][]) sq.getSqueeze().get(0); 

        sq.cal(MMCal.sum(inter2));
        double[] inter1 = (double[]) sq.getSqueeze().get(0);   
        
        int x0 = ((int[]) MMCal.maxWithIndex(inter1).get(1))[0];
        Assert.assertEquals(3, x0);
	}
	
	
	@Test
	public void y0Test(){	
    	double[][][] Hb = H.clone();   
    	double maxhi =  MMCal.max(MMCal.concatAll(H))[0];
    	Hb = MMCal.setValueWithFilter(0.0, H, MMCal.lt(0.5*maxhi,H)); 
		Squeeze sq = new MMCal.Squeeze();	
		sq.cal(MMCal.sum(Hb,1));
	    	    
	    double[][] inter2 = (double[][]) sq.getSqueeze().get(0);
        
	    sq.cal(MMCal.sum(MMCal.nonConjugateTranspose(inter2)));
    	double[] inter1 = (double[]) sq.getSqueeze().get(0);       

    	int y0 = ((int[]) MMCal.maxWithIndex(inter1).get(1))[0]; 
    	Assert.assertEquals(1, y0); 
	}
	
	@Test
	public void z0Test(){	
    	double[][][] Hb = H.clone();   
    	double maxhi =  MMCal.max(MMCal.concatAll(H))[0];
    	Hb = MMCal.setValueWithFilter(0.0, H, MMCal.lt(0.5*maxhi,H)); 
		Squeeze sq = new MMCal.Squeeze();	
	    sq.cal(MMCal.sum(Hb,0));
	    	    
	    double[][] inter2 = (double[][]) sq.getSqueeze().get(0);
	    sq.cal(MMCal.sum(inter2));    
    	double[] inter1 = (double[]) sq.getSqueeze().get(0);   
    	int z0 = ((int[]) MMCal.maxWithIndex(inter1).get(1))[0]; 
    	Assert.assertEquals(2, z0);
	}

	@Test
	public void rTest(){	
        int M2 = H[0].length;
        int M1 = H[0][0].length;
        int M3 = H.length; 
        int[] sizeOut = {5,3,3};
		int rX = Math.min(Math.floorDiv(M1,2),(sizeOut[0]-1)/2);
        int rY = Math.min(Math.floorDiv(M2,2),(sizeOut[1]-1)/2);  
        int rZ = Math.min(Math.floorDiv(M3,2),(sizeOut[2]-1)/2); 
        
    	Assert.assertEquals(2, rX);
    	Assert.assertEquals(1, rY);
    	Assert.assertEquals(1, rZ);
    	int x0 = 3;
        int[] r0X = new int[x0+rX-(x0-rX)+1];
        for (int i=0 ; i< r0X.length ; i++)
        	r0X[i] = x0-rX + i;
        int y0 = 1;
        int[] r0Y = new int[y0+rY-(y0-rY)+1];
        for (int i=0 ; i< r0Y.length ; i++)
        	r0Y[i] = y0-rY + i;
        int z0 = 2;
        int[] r0Z = new int[z0+rZ-(z0-rZ)+1];
        for (int i=0 ; i< r0Z.length ; i++)
        	r0Z[i] = z0-rZ + i;
        
        int[] expectedindX = {1,2,3,4,5};
        int[] expectedindY = {0,1,2};
        int[] expectedindZ = {1,2,3};    
        int[] indX = MMCal.max(MMCal.min(r0X,M1),0);       
        int[] indY = MMCal.max(MMCal.min(r0Y,M2),0);
        int[] indZ = MMCal.max(MMCal.min(r0Z,M3),0);
        Assert.assertArrayEquals(expectedindX, indX);
        Assert.assertArrayEquals(expectedindY, indY);
        Assert.assertArrayEquals(expectedindZ, indZ);
		
        // Remove negative values and normalize
        double[][][] Hnew = MMCal.getIndex(H, indX, indY, indZ);
        Hnew = MMCal.max(Hnew, 0.0);
        double normH = MMCal.sum( MMCal.concatAll(Hnew))[0];
        Hnew  = MMCal.rDivide ( Hnew , normH ); 
        double[][][] Hexpected = 
        	{{{   1.2308e-03 ,  6.1538e-04  , 1.2308e-03  , 1.8462e-03 ,  6.1538e-04},
        	  {   6.1538e-02 ,  1.2308e-01  , 1.8462e-01  , 6.1538e-02  , 1.2308e-03},
        	  {	   6.1538e-04 ,  1.2308e-02 ,  1.8462e-03,   6.1538e-04 ,  6.1538e-04}},
        	  {{   6.1538e-04 ,  1.2308e-03 ,  1.8462e-03 ,  6.1538e-04  , 6.1538e-04},
        	   {   6.1538e-02 ,  1.2308e-01  , 2.4615e-01  , 6.1538e-02 ,  1.2308e-03},
        	   {  6.1538e-04  , 1.2308e-02  , 1.8462e-03 ,  6.1538e-04  , 6.1538e-04}},
              {{ 6.1538e-04 ,  1.2308e-03 ,  1.8462e-03 ,  6.1538e-04 ,  6.1538e-04},
        	   {   6.1538e-03  , 1.2308e-03 ,  1.8462e-03 ,  6.1538e-04  , 1.2308e-03},
        	   {   6.1538e-04 ,  1.2308e-02 ,  1.8462e-03  , 6.1538e-04  , 6.1538e-04}}};
    	for (int k=0 ; k< Hnew.length ; k++)
    		for(int j=0 ; j< Hnew[0].length ; j++)
    			Assert.assertArrayEquals(Hexpected[k][j], Hnew[k][j],10E-6);

	}
	
	@Test
	public void testgetIndex(){
		double[][][] Hexpected 		
		 = {
        {{ 0.020, 0.01, 0.020, 0.03, 0.01},
        		{ 1.0, 2.0, 3.0, 1.0, 0.020},
        		{ 0.01, 0.2, 0.03, 0.01, 0.01}},
        {{ 0.01, 0.020, 0.03, 0.01, 0.01},
        			{  1.0, 2.0, 4.0, 1.0, 0.020},
        			{0.01, 0.2 ,0.03, 0.01, 0.01}},
        	{{ 0.01, 0.020, 0.03, 0.01, 0.01},
        			{ 0.1, 0.02, 0.03 ,0.01, 0.020},
        			{ 0.01, 0.2, 0.03, 0.01, 0.01}}};
		
		int[] index0 =  {1,2,3,4,5};
		int[] index1 = {0,1,2};
		int[] index2 = {1,2,3}; 
		double[][][] m = MMCal.getIndex(H, index0, index1, index2);
		Assert.assertArrayEquals(Hexpected, m);


		}
	
	@Test
	public void testTodoubleReal(){
		double[][][] mexpected = {{{1.0,0, 2.0,0, 3.0,0},{1.0,0, 2.0,0, 3.0,0},{1.0,0, 2.0,0, 3.0,0}},
		        {{4.0,0, 5.0,0, 6.0,0},{4.0,0, 5.0,0, 6.0,0},{4.0,0, 5.0,0, 6.0,0}},
		        {{5.0,0, 5.0,0, 5.0,0},{5.0,0, 5.0,0, 5.0,0},{5.0,0, 5.0,0, 5.0,0}}};
		double[][][] m = MMCal.complexCopy(mat333);
		Assert.assertArrayEquals(mexpected , m);
		
		double[][][] mm = MMCal.toDoubleReal(m); 
		Assert.assertArrayEquals(mat333 , mm);
	}
	
	@Test
	public void testconj(){
		double[][][] m = {{{1.0,-3.0, 2.0,2.0, 3.0,0},{1.0,0, 2.0,0, 3.0,0},{1.0,0, 2.0,0, 3.0,0}},
		        {{4.0,0, 5.0,0, 6.0,0},{4.0,6.2, 5.0,0, 6.0,0},{4.0,0, 5.0,0, 6.0,0}},
		        {{5.0,0, 5.0,0, 5.0,0},{5.0,-15.3, 5.0,0, 5.0,0},{5.0,0, 5.0,0.2, 5.0,0}}};
		m = MMCal.conj(m);
		double[][][] mexpected = {{{1.0,3.0, 2.0,-2.0, 3.0,0},{1.0,0, 2.0,0, 3.0,0},{1.0,0, 2.0,0, 3.0,0}},
		        {{4.0,0, 5.0,0, 6.0,0},{4.0,-6.2, 5.0,0, 6.0,0},{4.0,0, 5.0,0, 6.0,0}},
		        {{5.0,0, 5.0,0, 5.0,0},{5.0,15.3, 5.0,0, 5.0,0},{5.0,0, 5.0,-0.2, 5.0,0}}};
    	for (int k=0 ; k< m.length ; k++)
    		for(int j=0 ; j< m[0].length ; j++)
    			Assert.assertArrayEquals(mexpected[k][j] , m[k][j] ,10-3);
		
	}
	
	
	@Test
	public void HreductionTest(){
		double[][][] Hexpected = {{{1.2308e-03  , 6.1538e-04  , 1.2308e-03  , 1.8462e-03 ,  6.1538e-04},		
									{6.1538e-02  , 1.2308e-01,  1.8462e-01 ,  6.1538e-02  , 1.2308e-03},
									{ 6.1538e-04 ,  1.2308e-02 ,  1.8462e-03 ,  6.1538e-04 ,  6.1538e-04}},
								{{6.1538e-04 ,  1.2308e-03 ,  1.8462e-03  , 6.1538e-04   ,6.1538e-04},
									{6.1538e-02  , 1.2308e-01  ,2.4615e-01 ,  6.1538e-02  , 1.2308e-03},
									{6.1538e-04  , 1.2308e-02  , 1.8462e-03  , 6.1538e-04  , 6.1538e-04}},
								{{ 6.1538e-04 ,  1.2308e-03 ,  1.8462e-03  , 6.1538e-04 ,  6.1538e-04},
									{6.1538e-03 ,  1.2308e-03  , 1.8462e-03 ,  6.1538e-04  , 1.2308e-03 },
									 {6.1538e-04  , 1.2308e-02 ,  1.8462e-03 ,  6.1538e-04  , 6.1538e-04}}};
						
    	double[][][] Hb = H.clone();   
    	double maxhi =  MMCal.max(MMCal.concatAll(H))[0];
    	Hb = MMCal.setValueWithFilter(0.0, H, MMCal.lt(0.5*maxhi,H)); 
    	Assert.assertNotEquals(Hb, H);
		Squeeze sq = new MMCal.Squeeze();
		
		sq.cal(MMCal.sum(Hb,2));    
	    double[][] inter2 = (double[][]) sq.getSqueeze().get(0);
	    sq.cal(MMCal.sum(inter2));
    	double[] inter1 = (double[]) sq.getSqueeze().get(0);        
    	int y0 = ((int[]) MMCal.maxWithIndex(inter1).get(1))[0]; 

	    sq.cal(MMCal.sum(Hb,0));	    
	    inter2 = (double[][]) sq.getSqueeze().get(0);
	    sq.cal(MMCal.sum(inter2));
    	inter1 = (double[]) sq.getSqueeze().get(0);        
    	int z0 = ((int[]) MMCal.maxWithIndex(inter1).get(1))[0]; 
    	
        sq.cal(MMCal.sum(Hb,1));       
        inter2= (double[][]) sq.getSqueeze().get(0);       
        sq.cal(MMCal.sum(MMCal.nonConjugateTranspose(inter2)));
        inter1 = (double[]) sq.getSqueeze().get(0);   
        int x0 = ((int[]) MMCal.maxWithIndex(inter1).get(1))[0];
        Assert.assertEquals(1, x0);
        
		int M1 = H[0].length;
        int M2 = H[0][0].length;
        int M3 = H.length;        
        int[] sizeOut ={ 3,5,3};
        
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

        int[] indX = MMCal.max(MMCal.min(r0X,M1),0);       
        int[] indY = MMCal.max(MMCal.min(r0Y,M2),0);
        int[] indZ = MMCal.max(MMCal.min(r0Z,M3),0);
        Assert.assertArrayEquals(new int[]{0, 1,2},indX);
        Assert.assertArrayEquals(new int[]{1, 2,3,4,5},indY);
        Assert.assertArrayEquals(new int[]{1, 2,3},indZ);
        double[][][] Hnew = MMCal.getIndex(H, indY,indX ,indZ);
        // Remove negative values and normalize
        Hnew = MMCal.max(Hnew,0);
        
      //  System.out.print("sum Hnew "+Double.toString(MMCal.sum( MMCal.concatAll(Hnew))[0]));
        double normH = MMCal.sum( MMCal.concatAll(Hnew))[0];
        Hnew = MMCal.rDivide ( Hnew, normH ); 
        
        //New Matrix size
        M1 = Hnew[0].length;
        M2 = Hnew[0][0].length;
        M3 = Hnew.length;  
        Assert.assertEquals(Hexpected[0].length, Hnew[0].length);
        Assert.assertEquals(Hexpected[0][0].length, Hnew[0][0].length);
        Assert.assertEquals(Hexpected.length, Hnew.length);
    	for (int k=0 ; k< Hnew.length ; k++)
    		for(int j=0 ; j< Hnew[0].length ; j++)
    			Assert.assertArrayEquals(Hexpected[k][j], Hnew[k][j],10E-4);
	}
		
	private static String toString(int[] A){
		String t ="[";
		for(int i=0 ; i < A.length ; i++){ 
				t = t +  String.valueOf(A[i]) + " ";
			}
		t = t + "]";
		return t;
	}
	
	@Test
	public void interpCubic3Test(){
		double[][][] expected = 
		{{{ 0.010000  , 0.017214 ,  0.012875,   0.029193,   0.010000},
		   {0.010000 ,  0.043727  , 0.078969 ,  0.030804 ,  0.010000}}};
	    double[] resIN = {0.008, 0.008, 0.01};
	    double[] resOUT = {0.01, 0.01, 0.08};
    	double taillex = resIN[0]*H[0][0].length;
    	double tailley = resIN[1]*H[0].length;
    	double taillez = resIN[2]*H.length;
    	double[] x = new  double[H[0][0].length];
    	double[] y = new  double[H[0].length];
    	double[] z = new  double[H.length];
    	int taillexi = (int) (Math.floor(resIN[0]*(double)(H[0][0].length-1)/resOUT[0]))+1;
    	int tailleyi = (int) (Math.floor(resIN[1]*(double)(H[0].length-1)/resOUT[1]))+1;
    	int taillezi = (int) (Math.floor(resIN[2]*(double)(H.length-1)/resOUT[2]))+1;
    	double[] xi = new  double[taillexi];
    	double[] yi = new  double[tailleyi];
    	double[] zi = new  double[taillezi];
   	
    	for (int i=0; i < H[0][0].length; i++ ){
    		x[i] = i* resIN[0] ;    			
    	}
    	for (int i=0; i < H[0].length; i++ ){
    		y[i] = i* resIN[1] ;
    	}
    	for (int i=0; i < H.length; i++ ){
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
   //     System.out.print("x :"+toString(x));
   //     System.out.print("y :"+toString(y));
   //     System.out.print("z :"+toString(z));
   //     System.out.print("xi :"+toString(xi));
   //     System.out.print("yi :"+toString(yi));
   //     System.out.print("zi :"+toString(zi));
		double[][][] Vi = MMCal.interpCubic3(z,y,x,zi,yi,xi,H);
	//	System.out.print("Vi :"+toString(Vi));
		//for (int)
	//	for (int k=0 ; k< expected.length ; k++)
    //		for(int j=0 ; j< expected[0].length ; j++)
	//	Assert.assertArrayEquals(expected[k][j], Vi[k][j],0.01)
    			//;
		
	//	TricubicInterpolator tri = new TricubicInterpolator();
	//	TriCubicInterpolator triCub = new TriCubicInterpolator(z, y, x, H);
		
	//	TricubicInterpolatingFunction f = tri.interpolate(z, y, x, H);
	//	double v = f.value(zi[0], yi[0], xi[0]);
	//	Vi = triCub.interpolate(zi, yi, xi);
		for (int k=0 ; k< expected.length ; k++)
    		for(int j=0 ; j< expected[0].length ; j++){
	         Assert.assertArrayEquals(expected[k][j], Vi[k][j],5e-2);
	//         System.out.print("expected :"+toString(expected[k][j]));
	//         System.out.print("Vi :"+toString(Vi[k][j]));
    		}
	}
	
	@Test
	public void conv2Dadjoint_fourier(){
		double[][] mcal = MMCal.conv2Dadjoint_fourier(mat32, mat32);
		Assert.assertArrayEquals(mcal,mat32);
	}
	
	@Test
	public void testconv3D_fourier(){
	//	double[][][]  ones = new double[3][2][4];
	//	ones[0][1][2] =1.0;
    //    DoubleFFT_3D FFTCal1 = new DoubleFFT_3D(ones.length,ones[0].length,ones[0][0].length);
    //    double[][][] Dextcomplex = MMCal.complexCopy(ones);
     //   FFTCal1.complexForward(Dextcomplex);
     //   System.out.print("Dextcomplex fft :"+toString(Dextcomplex)); 
     //   DoubleFFT_3D FFTCal2 = new DoubleFFT_3D(ones.length,ones[0].length,ones[0][0].length);
    //    FFTCal2.complexInverse(Dextcomplex, true);
     //   System.out.print("Dextcomplex ifft :"+toString(Dextcomplex)); 
		double[][][] mcal = MMCal.conv3D_fourier(H, mat333);
		double[][][] mexpected = {{
			{  1.3200  ,  4.8700  , 11.0500  , 14.2000 ,  11.7700  ,  3.6100},
			{   1.4800 ,   5.3200 ,  12.4800 ,  16.0500,  13.2800 ,   3.8900},
			{   1.3100  ,  5.0400 ,  12.1300 ,  15.6700 , 12.8400  , 3.5500}},
		{{   5.7100,   18.4800 ,  40.3000 ,  47.4600 ,  38.1300  , 9.8200},
		{    6.1600  , 20.1100  , 43.0900 ,  50.7200  , 39.7200  ,10.1800},
		{    5.8500 ,  19.6400 ,  42.4900 ,  50.0100  , 38.9900  ,  9.7000}},
		{{  9.9400  , 29.3600 ,  63.0000 ,  66.8700  , 50.0000  , 11.6300},
		{   10.7900 ,  32.4300 ,  65.8400 ,  70.1300  , 50.6400 , 11.8900},
		{  10.5200  , 31.9700  , 65.2100  , 69.3900  , 49.9500  , 11.5300}},
		{{    6.1400 ,  16.5300  , 36.4000 ,  35.9200 ,  25.9300 ,   5.4700},
		{    7.2300 ,  19.5300 , 38.7800 ,  38.5100 ,  26.4500  ,  5.6800},
		{    7.0400 ,  19.1400 , 38.2000  ,37.9000 , 25.9300  ,  5.4700}}};
		
		for (int k=0 ; k< mcal.length ; k++)
    		for(int j=0 ; j< mcal[0].length ; j++)
    			Assert.assertArrayEquals(mexpected[k][j],mcal[k][j],10E-4);
	}
	
	
	@Test
	public void testFreqfilt3D(){
		int[] p = new int[3];

		double[][][] mm = new double[3][3][5];
		mm = MMCal.add(mm, 0.5);
		p[0] = Math.floorDiv(mm[0][0].length,2);
		p[1] = Math.floorDiv(mm[0].length,2);
		p[2] = Math.floorDiv(mm.length,2);
		/*ComplexNum[][][] ccexpected =  
			{{{	new DoublePrecComplNum( 22.50000 , 0.00000),  new DoublePrecComplNum(12.95723 , 0.00000), new DoublePrecComplNum(-2.39440 ,0.00000),new DoublePrecComplNum(-4.50000,0.00000),
					new DoublePrecComplNum( 2.93717 ,-  0.00000), new DoublePrecComplNum( 2.93717 ,-  0.00000),  new DoublePrecComplNum(-4.50000 ,  0.00000),new DoublePrecComplNum(-2.39440 ,-  0.00000),new DoublePrecComplNum( 12.95723 ,  0.00000)},												                      
			  { new DoublePrecComplNum(12.13525 , 0.00000),  new DoublePrecComplNum( 6.98841 ,0.00000), new DoublePrecComplNum(-1.29141 , 0.00000),new DoublePrecComplNum(-2.42705 ,  0.00000),
					new DoublePrecComplNum(1.58415 ,-  0.00000),new DoublePrecComplNum(  1.58415 ,  0.00000),new DoublePrecComplNum(-2.42705 , 0.00000),new DoublePrecComplNum(-1.29141 ,  0.00000),new DoublePrecComplNum( 6.98841 ,-  0.00000)},
			  {	new DoublePrecComplNum(-4.63525 ,  0.00000), new DoublePrecComplNum(-2.66934 ,  0.00000),new DoublePrecComplNum(0.49327,  0.00000),new DoublePrecComplNum(  0.92705 ,  0.00000),
					new DoublePrecComplNum(1.58415 ,-  0.00000),new DoublePrecComplNum(-0.60509 ,  0.00000),new DoublePrecComplNum( 0.92705 ,  0.00000),new DoublePrecComplNum( 0.49327 ,-  0.00000),new DoublePrecComplNum(-2.66934 ,-  0.00000)},
			  {	new DoublePrecComplNum(-4.63525 ,  0.00000), new DoublePrecComplNum(-2.66934 ,  0.00000),new DoublePrecComplNum(0.49327 ,  0.00000),new DoublePrecComplNum( 0.92705 ,-  0.00000),
					new DoublePrecComplNum(-0.60509 ,-  0.00000),new DoublePrecComplNum(-0.60509 ,  0.00000),new DoublePrecComplNum( 0.92705 ,-  0.00000),new DoublePrecComplNum( 0.49327 ,-  0.00000),new DoublePrecComplNum(-2.66934 ,-  0.00000)},
			  {	new DoublePrecComplNum(12.13525 ,  0.00000), new DoublePrecComplNum( 6.98841 ,  0.00000),new DoublePrecComplNum( -1.29141 ,-  0.00000),new DoublePrecComplNum(-2.42705 ,-  0.00000),
					new DoublePrecComplNum( 1.58415 ,-  0.00000),new DoublePrecComplNum( 1.58415 ,  0.00000),new DoublePrecComplNum(-2.42705 ,-  0.00000),new DoublePrecComplNum(-1.29141, -  0.00000),new DoublePrecComplNum( 6.98841 ,-  0.00000)}},					  				    		  
			  
			 {{	new DoublePrecComplNum(12.13525 ,  0.00000),  new DoublePrecComplNum( 6.98841 ,  0.00000),new DoublePrecComplNum(-1.29141 ,  0.00000),new DoublePrecComplNum( -2.42705 +  0.00000),
				    new DoublePrecComplNum(1.58415 ,  0.00000),new DoublePrecComplNum( 1.58415 ,-  0.00000),new DoublePrecComplNum( -2.42705 ,  0.00000),new DoublePrecComplNum(-1.29141 ,-  0.00000),new DoublePrecComplNum(6.98841 ,  0.00000)},
			  {	new DoublePrecComplNum(6.54508 ,  0.00000),  new DoublePrecComplNum( 3.76916 , 0.00000),new DoublePrecComplNum(-0.69651 ,  0.00000),new DoublePrecComplNum(-1.30902 +  0.00000),
				    new DoublePrecComplNum( 0.85440 ,-  0.00000),new DoublePrecComplNum(0.85440 ,  0.00000),new DoublePrecComplNum( -1.30902 ,  0.00000),new DoublePrecComplNum(-0.69651 ,  0.00000),new DoublePrecComplNum(3.76916 ,-  0.00000)},
			  {	new DoublePrecComplNum(-2.50000 ,  0.00000),  new DoublePrecComplNum(-1.43969 ,  0.00000),new DoublePrecComplNum( 0.26604 ,  0.00000),new DoublePrecComplNum( 0.50000 +  0.00000),new DoublePrecComplNum( -0.32635 ,-  0.00000),new DoublePrecComplNum(-0.32635 ,  0.00000),new DoublePrecComplNum(0.50000 ,  0.00000),new DoublePrecComplNum(0.26604 ,-  0.00000),new DoublePrecComplNum(-1.43969 ,-  0.00000)},
			  {	new DoublePrecComplNum(-2.50000 , -  0.00000),new DoublePrecComplNum(-1.43969 , 0.00000),
					new DoublePrecComplNum(0.26604 ,  0.00000),new DoublePrecComplNum(  0.50000 -  0.00000),new DoublePrecComplNum(-0.32635 ,-  0.00000),new DoublePrecComplNum( -0.32635 ,  0.00000),new DoublePrecComplNum(0.50000 ,-  0.00000),new DoublePrecComplNum(0.26604 ,-  0.00000),new DoublePrecComplNum(-1.43969 ,-  0.00000)},
			  {	new DoublePrecComplNum(6.54508 ,-  0.00000),  new DoublePrecComplNum(3.76916 ,  0.00000),
					new DoublePrecComplNum( -0.69651 ,-  0.00000),new DoublePrecComplNum(-1.30902 -  0.00000),new DoublePrecComplNum( 0.85440 ,-  0.00000),new DoublePrecComplNum(0.85440 ,  0.00000),new DoublePrecComplNum(-1.30902 ,-  0.00000),new DoublePrecComplNum(-0.69651 ,-  0.00000),new DoublePrecComplNum(3.76916 ,-  0.00000)}},
				   
			  
			 {{ new DoublePrecComplNum(-4.63525 , 0.00000), new DoublePrecComplNum(-2.66934 , 0.00000),new DoublePrecComplNum(0.49327 ,- 0.00000),new DoublePrecComplNum( 0.92705 + 0.00000),new DoublePrecComplNum(-0.60509 ,- 0.00000),new DoublePrecComplNum(-0.60509 , 0.00000),new DoublePrecComplNum(0.92705 + 0.00000),new DoublePrecComplNum(  0.49327 + 0.00000),new DoublePrecComplNum( -2.66934 + 0.00000)},
			  {	new DoublePrecComplNum(-2.50000 , 0.00000),new DoublePrecComplNum(-1.43969 ,- 0.00000),new DoublePrecComplNum(0.26604 , 0.00000),new DoublePrecComplNum( 0.50000 + 0.00000),new DoublePrecComplNum( -0.32635 , 0.00000),new DoublePrecComplNum( -0.32635 ,- 0.00000),new DoublePrecComplNum( 0.50000 + 0.00000),new DoublePrecComplNum(   0.26604 + 0.00000),new DoublePrecComplNum(  -1.43969 + 0.00000)},
			  {	new DoublePrecComplNum(0.95492 , 0.00000), new DoublePrecComplNum(0.54991 ,- 0.00000), new DoublePrecComplNum(-0.10162 ,- 0.00000),new DoublePrecComplNum(-0.19098 + 0.00000),new DoublePrecComplNum( 0.12466 , 0.00000),new DoublePrecComplNum(0.12466 ,- 0.00000),new DoublePrecComplNum( -0.19098 + 0.00000),new DoublePrecComplNum(  -0.10162 + 0.00000),new DoublePrecComplNum(   0.54991 + 0.00000)},
			  {	new DoublePrecComplNum(0.95492 ,- 0.00000), new DoublePrecComplNum(0.54991 ,- 0.00000),new DoublePrecComplNum(-0.10162 ,- 0.00000),new DoublePrecComplNum( -0.19098 - 0.00000),new DoublePrecComplNum( 0.12466 , 0.00000),new DoublePrecComplNum(0.12466 ,- 0.00000),new DoublePrecComplNum( -0.19098 - 0.00000),new DoublePrecComplNum(  -0.10162 + 0.00000),new DoublePrecComplNum(   0.54991 + 0.00000)},
			  {	new DoublePrecComplNum(-2.50000 ,- 0.00000),new DoublePrecComplNum(-1.43969, - 0.00000),new DoublePrecComplNum( 0.26604 ,- 0.00000),new DoublePrecComplNum( 0.50000 - 0.00000),new DoublePrecComplNum( -0.32635 , 0.00000),new DoublePrecComplNum( -0.32635 ,- 0.00000),new DoublePrecComplNum(0.50000 - 0.00000),new DoublePrecComplNum(  0.26604 - 0.00000),new DoublePrecComplNum(  -1.43969 + 0.00000)}},
					    
			 {{	new DoublePrecComplNum( -4.63525 , 0.00000),new DoublePrecComplNum( -2.66934 , 0.00000),new DoublePrecComplNum( 0.49327 ,- 0.00000),new DoublePrecComplNum( 0.92705 + 0.00000),
				    new DoublePrecComplNum(-0.60509 ,- 0.00000),new DoublePrecComplNum(-0.60509 , 0.00000),new DoublePrecComplNum(   0.92705 + 0.00000),new DoublePrecComplNum(   0.49327 + 0.00000),new DoublePrecComplNum(  -2.66934 + 0.00000)},
			  {	new DoublePrecComplNum( -2.50000 , 0.00000), new DoublePrecComplNum( -1.43969 ,- 0.00000), new DoublePrecComplNum( 0.26604 , 0.00000),new DoublePrecComplNum( 0.50000 + 0.00000),
					new DoublePrecComplNum(-0.32635 , 0.00000),new DoublePrecComplNum(-0.32635 ,- 0.00000),new DoublePrecComplNum(   0.50000 + 0.00000),new DoublePrecComplNum(   0.26604 + 0.00000),new DoublePrecComplNum( -1.43969 + 0.00000)},
			  {	new DoublePrecComplNum(	0.95492 , 0.00000),new DoublePrecComplNum( 0.54991 ,- 0.00000), new DoublePrecComplNum( -0.10162 ,- 0.00000),new DoublePrecComplNum( -0.19098 + 0.00000),
					new DoublePrecComplNum( 0.12466 , 0.00000),new DoublePrecComplNum(0.12466 ,- 0.00000),new DoublePrecComplNum(   -0.19098 + 0.00000),new DoublePrecComplNum(  -0.10162 + 0.00000),new DoublePrecComplNum(   0.54991 + 0.00000)},
			  {	new DoublePrecComplNum(	0.95492 ,- 0.00000), new DoublePrecComplNum( 0.54991 ,- 0.00000),  new DoublePrecComplNum( -0.10162 ,- 0.00000),new DoublePrecComplNum( -0.19098 - 0.00000),
					new DoublePrecComplNum(0.12466 , 0.00000),new DoublePrecComplNum(0.12466 ,- 0.00000),new DoublePrecComplNum(   -0.19098 - 0.00000),new DoublePrecComplNum( -0.10162 + 0.00000),new DoublePrecComplNum(   0.54991 + 0.00000)},
			  {	new DoublePrecComplNum(	-2.50000 ,- 0.00000),new DoublePrecComplNum( -1.43969 ,- 0.00000), new DoublePrecComplNum( 0.26604 ,- 0.00000),new DoublePrecComplNum( 0.50000 - 0.00000),
					new DoublePrecComplNum(-0.32635 , 0.00000),new DoublePrecComplNum(-0.32635 ,- 0.00000),new DoublePrecComplNum(    0.50000 - 0.00000),new DoublePrecComplNum(   0.26604 - 0.00000),new DoublePrecComplNum(  -1.43969 + 0.00000)}},
	 			    
			 {{	new DoublePrecComplNum(	12.13525 , 0.00000),new DoublePrecComplNum(  6.98841 ,  0.00000),  new DoublePrecComplNum(-1.29141 , 0.00000), new DoublePrecComplNum( -2.42705 ,  0.00000),
				    new DoublePrecComplNum(1.58415 ,  0.00000),new DoublePrecComplNum(1.58415 ,-  0.00000),new DoublePrecComplNum(-2.42705 ,  0.00000),new DoublePrecComplNum(-1.29141 ,-  0.00000),new DoublePrecComplNum( 6.98841 ,  0.00000)},
			  {	new DoublePrecComplNum(	6.54508 , 0.00000), new DoublePrecComplNum( 3.76916 ,  0.00000),   new DoublePrecComplNum(-0.69651 ,  0.00000), new DoublePrecComplNum(-1.30902 ,  0.00000),
				    new DoublePrecComplNum( 0.85440 ,-  0.00000),new DoublePrecComplNum(0.85440 ,  0.00000),new DoublePrecComplNum(-1.30902 ,  0.00000),new DoublePrecComplNum(-0.69651 ,  0.00000),new DoublePrecComplNum(3.76916 ,-  0.00000)},
			  {	new DoublePrecComplNum(	-2.50000 ,  0.00000),new DoublePrecComplNum( -1.43969 , 0.00000), new DoublePrecComplNum( 0.26604 ,  0.00000), new DoublePrecComplNum( 0.50000,  0.00000) ,
				    new DoublePrecComplNum(-0.32635 ,-  0.00000),new DoublePrecComplNum(-0.32635 ,  0.00000),new DoublePrecComplNum( 0.50000 ,  0.00000),new DoublePrecComplNum(0.26604 ,-  0.00000),new DoublePrecComplNum(-1.43969 ,-  0.00000)},
			  {	new DoublePrecComplNum(-2.50000 ,-  0.00000), new DoublePrecComplNum( -1.43969 ,  0.00000), new DoublePrecComplNum( 0.26604 ,  0.00000), new DoublePrecComplNum( 0.50000 ,-  0.00000),
				    new DoublePrecComplNum( -0.32635 ,- 0.00000),new DoublePrecComplNum(-0.32635 , 0.00000),new DoublePrecComplNum(0.50000 ,-  0.00000),new DoublePrecComplNum(0.26604 ,-  0.00000),new DoublePrecComplNum(-1.43969 ,-  0.00000)},
			  {	new DoublePrecComplNum(	6.54508 ,-  0.00000),new DoublePrecComplNum( 3.76916 , 0.00000),  new DoublePrecComplNum( -0.69651 ,-  0.00000),  new DoublePrecComplNum( -1.30902 ,-  0.00000),
				    new DoublePrecComplNum( 0.85440, -  0.00000),new DoublePrecComplNum( 0.85440 ,  0.00000),new DoublePrecComplNum(-1.30902 ,-  0.00000),new DoublePrecComplNum(-0.69651 ,-  0.00000),new DoublePrecComplNum(3.76916 -  0.00000)}}};*/
		double[][][] ccexpected =  
			{{{	22.50000 , 0.00000,   12.95723 , 0.00000, -2.39440 ,0.00000,-4.50000,0.00000,
					 2.93717 ,-0.00000,  2.93717 ,-  0.00000,  -4.50000 ,  0.00000,-2.39440 ,-  0.00000, 12.95723 ,  0.00000},												                      
			  { 12.13525 , 0.00000,   6.98841 ,0.00000, -1.29141 , 0.00000,-2.42705 ,  0.00000,
					1.58415 , - 0.00000,  1.58415 ,  0.00000,-2.42705 , 0.00000,-1.29141 ,  0.00000, 6.98841 ,-  0.00000},
			  {	-4.63525 ,  0.00000, -2.66934 ,  0.00000,0.49327,  0.00000,  0.92705 ,  0.00000,
					1.58415 ,-  0.00000,-0.60509 ,  0.00000, 0.92705 ,  0.00000, 0.49327 ,-  0.00000,-2.66934 ,-  0.00000},
			  {	-4.63525 ,  0.00000, -2.66934 ,  0.00000,0.49327 ,  0.00000, 0.92705 ,-  0.00000,
					-0.60509 ,-  0.00000,-0.60509 ,  0.00000, 0.92705 ,-  0.00000, 0.49327 ,-  0.00000,-2.66934 ,-  0.00000},
			  {	12.13525 ,  0.00000,  6.98841 ,  0.00000, -1.29141 ,-  0.00000,-2.42705 ,-  0.00000,
					 1.58415 ,-  0.00000, 1.58415 ,  0.00000,-2.42705 ,-  0.00000,-1.29141, -  0.00000, 6.98841 ,-  0.00000}},					  				    		  
			  
			 {{	12.13525 ,  0.00000,   6.98841 ,  0.00000,-1.29141 ,  0.00000, -2.42705 ,  0.00000,
				    1.58415 ,  0.00000, 1.58415 ,-  0.00000, -2.42705 ,  0.00000,-1.29141 ,-  0.00000,6.98841 ,  0.00000},
			  {	6.54508 ,  0.00000,   3.76916 , 0.00000,-0.69651 ,  0.00000,-1.30902 , 0.00000,
				     0.85440 ,-  0.00000,0.85440 ,  0.00000, -1.30902 ,  0.00000,-0.69651 ,  0.00000,3.76916 ,-  0.00000},
			  {	-2.50000 ,  0.00000,  -1.43969 ,  0.00000, 0.26604 ,  0.00000, 0.50000 ,  0.00000, -0.32635 ,-  0.00000,-0.32635 ,  0.00000,0.50000 ,  0.00000,0.26604 ,-  0.00000,-1.43969 ,-  0.00000},
			  {	-2.50000 , -  0.00000,-1.43969 , 0.00000,
					0.26604 ,  0.00000,  0.50000 ,-  0.00000,-0.32635 ,-  0.00000, -0.32635 ,  0.00000,0.50000 ,-  0.00000,0.26604 ,-  0.00000,-1.43969 ,-  0.00000},
			  {	6.54508 ,-  0.00000,  3.76916 ,  0.00000,
					 -0.69651 ,-  0.00000,-1.30902, -  0.00000, 0.85440 ,-  0.00000,0.85440 ,  0.00000,-1.30902 ,-  0.00000,-0.69651 ,-  0.00000,3.76916 ,-  0.00000}},
				   
			  
			 {{ -4.63525 , 0.00000, -2.66934 , 0.00000,0.49327 ,- 0.00000, 0.92705 , 0.00000,-0.60509 ,- 0.00000,-0.60509 , 0.00000,0.92705 , 0.00000,  0.49327 , 0.00000, -2.66934, 0.00000},
			  {	-2.50000 , 0.00000,-1.43969 ,- 0.00000,0.26604 , 0.00000, 0.50000 , 0.00000, -0.32635 , 0.00000, -0.32635 ,- 0.00000, 0.50000 , 0.00000,   0.26604 , 0.00000,  -1.43969 , 0.00000},
			  {	0.95492 , 0.00000, 0.54991 ,- 0.00000, -0.10162 ,- 0.00000,-0.19098 , 0.00000, 0.12466 , 0.00000,0.12466 ,- 0.00000, -0.19098 , 0.00000,  -0.10162 , 0.00000,   0.54991 , 0.00000},
			  {	0.95492 ,- 0.00000, 0.54991 ,- 0.00000,-0.10162 ,- 0.00000, -0.19098, - 0.00000, 0.12466 , 0.00000,0.12466 ,- 0.00000, -0.19098 ,- 0.00000,  -0.10162 , 0.00000,   0.54991 , 0.00000},
			  {	-2.50000 ,- 0.00000,-1.43969, - 0.00000, 0.26604 ,- 0.00000, 0.50000 ,- 0.00000, -0.32635 , 0.00000, -0.32635 ,- 0.00000,0.50000, - 0.00000,  0.26604 ,- 0.00000,  -1.43969 , 0.00000}},
					    
			 {{	 -4.63525 , 0.00000, -2.66934 , 0.00000, 0.49327 ,- 0.00000, 0.92705 , 0.00000,
				    -0.60509 ,- 0.00000,-0.60509 , 0.00000,   0.92705 , 0.00000,   0.49327 , 0.00000,  -2.66934 , 0.00000},
			  {	 -2.50000 , 0.00000,  -1.43969 ,- 0.00000,  0.26604 , 0.00000, 0.50000 , 0.00000,
					-0.32635 , 0.00000,-0.32635 ,- 0.00000,   0.50000 , 0.00000,   0.26604 , 0.00000, -1.43969 , 0.00000},
			  {		0.95492 , 0.00000, 0.54991 ,- 0.00000,  -0.10162 ,- 0.00000, -0.19098 , 0.00000,
					 0.12466 , 0.00000,0.12466 ,- 0.00000,   -0.19098 , 0.00000,  -0.10162 , 0.00000,   0.54991 , 0.00000},
			  {		0.95492 ,- 0.00000,  0.54991 ,- 0.00000,   -0.10162 ,- 0.00000, -0.19098, - 0.00000,
					0.12466 , 0.00000,0.12466 ,- 0.00000,   -0.19098 ,- 0.00000, -0.10162 , 0.00000,   0.54991 , 0.00000},
			  {		-2.50000 ,- 0.00000, -1.43969 ,- 0.00000,  0.26604 ,- 0.00000, 0.50000, - 0.00000,
					-0.32635 , 0.00000,-0.32635 ,- 0.00000,    0.50000 ,- 0.00000,   0.26604 ,- 0.00000,  -1.43969 , 0.00000}},
	 			    
			 {{		12.13525 , 0.00000,  6.98841 ,  0.00000,  -1.29141 , 0.00000,  -2.42705 ,  0.00000,
				    1.58415 ,  0.00000,1.58415 ,-  0.00000,-2.42705 ,  0.00000,-1.29141 ,-  0.00000, 6.98841 ,  0.00000},
			  {		6.54508 , 0.00000,  3.76916 ,  0.00000,   -0.69651 ,  0.00000, -1.30902 ,  0.00000,
				     0.85440 ,-  0.00000,0.85440 ,  0.00000,-1.30902 ,  0.00000,-0.69651 ,  0.00000, 3.76916 ,-  0.00000},
			  {		-2.50000 ,  0.00000, -1.43969 , 0.00000,  0.26604 ,  0.00000,  0.50000,  0.00000 ,
				    -0.32635 ,-  0.00000,-0.32635 ,  0.00000, 0.50000 ,  0.00000,0.26604 ,-  0.00000,-1.43969 ,-  0.00000},
			  {	-2.50000 ,-  0.00000,  -1.43969 ,  0.00000,  0.26604 ,  0.00000,  0.50000 ,-  0.00000,
				     -0.32635 ,- 0.00000,-0.32635 , 0.00000,0.50000 ,-  0.00000,0.26604 ,-  0.00000,-1.43969 ,-  0.00000},
			  {		6.54508 ,-  0.00000, 3.76916 , 0.00000,   -0.69651 ,-  0.00000,   -1.30902 ,-  0.00000,
				     0.85440, -  0.00000, 0.85440 ,  0.00000,-1.30902 ,-  0.00000,-0.69651 ,-  0.00000,3.76916 ,-  0.00000}}};		    
	//	System.out.print("n "+ String.valueOf(2*p[0]+mm[0][0].length)  );
	//	System.out.print("m "+ String.valueOf(2*p[1]+mm[0].length)  );
	//	System.out.print("p "+ String.valueOf(2*p[2]+mm.length)  );
		double[][][] cc = MMCal.freqfilt3D(mm,
				2*p[0]+mm[0][0].length, 2*p[1]+mm[0].length,2*p[2]+mm.length );
    	for (int k=0 ; k< cc.length ; k++)
    		for(int j=0 ; j< cc[0].length ; j++)
    			for(int i=0 ; i< cc[0][0].length/2 ; i++){
    				
    				if (!(k==0 && (j == 4 || j==2)&& i==4))
    				{
    					// System.out.print(" i: "+String.valueOf(i));
    					Assert.assertEquals(ccexpected[k][j][2*i],cc[k][j][2*i],10E-1);
    				}
    			//	System.out.print(String.valueOf(ccexpected[k][j][2*i+1]));
    			//	System.out.print(String.valueOf(cc[k][j][2*i+1]));
    				Assert.assertEquals(ccexpected[k][j][2*i+1],cc[k][j][2*i+1],10E-4);
    			}
	}
	
	private static String toString(double[] A){
		String t ="[";
			for(int i=0 ; i < A.length ; i++){ 
				t = t +  String.valueOf(A[i]) + " ";
			}
		t = t + "]";
		return t;
	}
	private static String toString(double[][] A){
		String t ="[";
		for(int j=0 ; j < A.length ; j++)
		{
			t = t + "[";
			for(int i=0 ; i < A[0].length ; i++){ 
				t = t +  String.valueOf(A[j][i]) + " ";
			}
			t = t + "]";
		}
		t = t + "]";
		return t;
	}
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
				t = t + "]";
		}
			t = t + "]";
		}
		t = t + "]";
		return t;
	}	
	//public static void main(String[] args) {
//		org.junit.runner.JUnitCore.main(MMCalTest);
//		}

}
