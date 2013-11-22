inp="/Volumes/Temp/NicolaDM/PoMo_Phylogenetic_project/HyPhy_concatenation/ILS_10Ne_10S_20G_R1.txt";
out2="/Volumes/Temp/NicolaDM/PoMo_Phylogenetic_project/HyPhy_concatenation/ILS_10Ne_10S_20G_R1_swap_HKY_out.txt";
DataSet ds = ReadDataFile(inp);
DataSetFilter filteredData = CreateFilter(ds,1,"","");
HarvestFrequencies (F, filteredData, 1, 1, 1);


/*Define global parameters*/
global kappa=1.0; 
global mu=1.0;

/*Define ancestral frequencies*/
Freqs={{F[0] }{ F[1] }{ F[2] }{ F[3] }};


/*Define matrix:*/ 
matrix1 = {{*, t * mu, t * kappa, t * mu}{t * mu, *, t * mu, t * kappa}{t * kappa, t * mu, *, t * mu}{t * mu, t * kappa, t * mu , *}};

Model M1 = (matrix1, Freqs, 1);
ACCEPT_ROOTED_TREES=1;
AUTOMATICALLY_CONVERT_BRANCH_LENGTHS = 1;

A_DISTANCE_METHOD 	   = 1;
SHORT_MPI_RETURN 	   = 1;
swapType=0;
dataType=0;
MESSAGE_LOGGING = 0;
totalRearrangementsTried = 0;
totalSupportWeight		 = 1;
_DO_TREE_REBALANCE_ = 1;
_KEEP_I_LABELS_ = 1;





/*NJ*/

distanceChoice=1;
A_DISTANCE_METHOD   = 1;
ACCEPT_ROOTED_TREES=1;
AUTOMATICALLY_CONVERT_BRANCH_LENGTHS = 1;
methodIndex=1;

DISTANCE_PROMPTS  = 1;
p = InferTreeTopology(1.0);
DISTANCE_PROMPTS = 0;

/* now with the treeNodes matrix ready we can convert it into a Newick string */

treeString = TreeMatrix2TreeString (0);/*1*/
fprintf (stdout, "\n\n --------------------- INFERRED TREE --------------------- \n\n", treeString,"\n");
fprintf (out2,CLEAR_FILE,treeString,";\n\n");













/*Find Root*/

Model M1 = (matrix1, Freqs, 0);
ACCEPT_ROOTED_TREES=1;
AUTOMATICALLY_CONVERT_BRANCH_LENGTHS = 1;

treeString = RerootTree (treeString, 0);

Tree   testTree2 = treeString;
Lea2= TipCount(testTree2);
Int2= BranchCount(testTree2);
BranchCount2= Lea2+Int2;

for (ibc2 = 0; ibc2 < BranchCount2; ibc2=ibc2+1)				/*iterating over roots*/
{

	Model M1 = (matrix1, Freqs, 0);
	if (ibc2<Int2)
	{
		Tree Tree2 = treeString;
		aBranch2   = BranchName (Tree2,ibc2);
		testTreest = RerootTree (Tree2, aBranch2);
		Tree testTree=testTreest;
		fprintf (stdout, "\n Internal node, tree:	",testTree,"\n");
	}
	else
	{
		Tree Tree2 = treeString;
		aBranch2   = TipName (Tree2,ibc2-Int2);
		testTreest = RerootTree (Tree2, aBranch2);
		Tree testTree=testTreest;
		fprintf (stdout, "\n External node, tree:	",testTree,"\n");
	}

	LikelihoodFunction theLF = (filteredData, testTree);
	MolecularClock (testTree,t);
	Optimize (res, theLF);
	fprintf (stdout, "\n One root, after optimization, lk:	",res[1][0],"\n");
	if (ibc2==0)
	{
		bestabs=res[1][0];
		fprintf (out2,theLF,"\n\n");
		fprintf (stdout, theLF,"\n\n");
		bestTree=testTreest;
	}

	if (res[1][0]>bestabs)
	{
		bestabs=res[1][0];
		fprintf (out2,theLF,"\n\n");
		fprintf (stdout, theLF,"\n\n");
		bestTree=testTreest;
	}

}




























function TreeMatrix2TreeString (doLengths)

{
	treeString = "";
	p = 0;
	k = 0;
	m = treeNodes[0][1];
	n = treeNodes[0][0];
	treeString*(Rows(treeNodes)*25);
	while (m)
	{	
		
		if (m>p)
		{
			if (p)
			{
				treeString*",";
			}
			for (j=p;j<m;j=j+1)
			{
				treeString*"(";
			}
		}
		else
		{
			if (m<p)
			{
				for (j=m;j<p;j=j+1)
				{
					treeString*")";
				}
			}	 
			else
			{
				treeString*",";
			}	 
		}
		
		if (n<ds.species)
		{
			GetString (nodeName, ds, n);
			if (doLengths != 1)
			{
				
				treeString*nodeName;	
						
			}
			else
			{	 
				
				treeString*taxonNameMap[nodeName];
				
			}
		}
		if (doLengths>.5)
		{
			nodeName = ":"+treeNodes[k][2];
			treeString*nodeName;
		}
		k=k+1;
		p=m;
		n=treeNodes[k][0];
		m=treeNodes[k][1];
	}
	for (j=m;j<p;j=j+1)
	{
		treeString*")";
	}
	treeString*0;
	return treeString;
}




























function InferTreeTopology(verbFlag)
{
	distanceMatrix = {ds.species,ds.species};

	if (distanceChoice)
	{
		if (distanceChoice == 2)
		{
			SetDialogPrompt ("Load the distance matrix");
			fscanf (PROMPT_FOR_FILE,"NMatrix",distanceMatrix);
			if ((Rows(distanceMatrix) != ds.species)||(Columns(distanceMatrix) != ds.species))
			{
				fprintf (stdout, "\nThe dimensions of the distance matrix are incompatible with the data set.\n");
				return  0;
			}
		}
		else
		{
			if (verbFlag)
			{
				fprintf (stdout, "\nHYPHY Kernel is computing pairwise maximum likelihood distance estimates. A total of ", Format(ds.species*(ds.species-1)/2,0,0), " estimations will be performed.\n");
			}

			_pddVF = 1;
			ExecuteAFile ("pairwiseDistanceEstimator.ibf");
		}
	}
	else
	{
		/*#include "chooseDistanceFormula.def";*/
		
		dummy = InitializeDistances (0);
		
		for (i = 0; i<ds.species; i=i+1)
		{
			for (j = i+1; j<ds.species; j = j+1)
			{
				distanceMatrix[i][j] = ComputeDistanceFormula (i,j);
			}
		}
	}

	MESSAGE_LOGGING 		 	= 1;
	cladesMade 					= 1;
	

	if (ds.species == 2)
	{
		d1 = distanceMatrix[0][1]/2;
		treeNodes = {{0,1,d1__},
					 {1,1,d1__},
					 {2,0,0}};
					 
		cladesInfo = {{2,0}};
	}
	else
	{
		if (ds.species == 3)
		{
			d1 = (distanceMatrix[0][1]+distanceMatrix[0][2]-distanceMatrix[1][2])/2;
			d2 = (distanceMatrix[0][1]-distanceMatrix[0][2]+distanceMatrix[1][2])/2;
			d3 = (distanceMatrix[1][2]+distanceMatrix[0][2]-distanceMatrix[0][1])/2;
			treeNodes = {{0,1,d1__},
						 {1,1,d2__},
						 {2,1,d3__}
						 {3,0,0}};
						 
			cladesInfo = {{3,0}};		
		}
		else
		{	
			njm = (distanceMatrix > methodIndex)>=ds.species;
				
			treeNodes 		= {2*(ds.species+1),3};
			cladesInfo	    = {ds.species-1,2};
			
			for (i=Rows(treeNodes)-1; i>=0; i=i-1)
			{
				treeNodes[i][0] = njm[i][0];
				treeNodes[i][1] = njm[i][1];
				treeNodes[i][2] = njm[i][2];
			}

			for (i=Rows(cladesInfo)-1; i>=0; i=i-1)
			{
				cladesInfo[i][0] = njm[i][3];
				cladesInfo[i][1] = njm[i][4];
			}
			
			njm = 0;
		}
	}

	fprintf (stdout, "\n\n --------------------- INFERRED MATRIX --------------------- \n\n", distanceMatrix);

	fprintf (stdout, "\n\n --------------------- INFERRED treeNodes --------------------- \n\n", treeNodes);

	distanceMatrix = 0;
	
	return 1.0;
}
