inp =  "/Volumes/Temp/NicolaDM/PoMo_Phylogenetic_project/HyPhy_concatenation/popul_10Ne_10S_200G_R1_PoMo9.txt";
out2 = "/Volumes/Temp/NicolaDM/PoMo_Phylogenetic_project/HyPhy_concatenation/popul_10Ne_10S_200G_R1_Nuc_NNIwithRoot_out.txt";
treeString="";
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













/*Preliminary ML*/


Tree givenTree=treeString;
LikelihoodFunction theLF = (filteredData, givenTree);
MolecularClock (givenTree,t);
Optimize (res, theLF);
fprintf (stdout, "\n\n --------------------- Preliminary ML --------------------- \n\n", theLF,"\n");
fprintf (out2,CLEAR_FILE,"Preliminary ML:\n",theLF,"\n\n");
BestTreeEver = treeString;
bestabs=res[1][0];
valueToBeat=res[1][0];







/* SWAPPING JOB*/

fprintf (stdout, "\n\n --------------------- START SWAPPING --------------------- \n\n");

A_DISTANCE_METHOD 	   = 1;
SHORT_MPI_RETURN 	   = 1;
swapType=0;
dataType=0;
#include "heuristicMethodNPBootstrap.bf";
MESSAGE_LOGGING = 0;
totalRearrangementsTried = 0;
totalSupportWeight		 = 1;
_DO_TREE_REBALANCE_ = 1;
_KEEP_I_LABELS_ = 1;

if (_KEEP_I_LABELS_)
{
	INTERNAL_NODE_PREFIX = "intNode";
}
	
Tree givenTree = treeString;
if (_KEEP_I_LABELS_)
{
	INTERNAL_NODE_PREFIX = "Node";
}

_KEEP_I_LABELS_ = 0;
globalOption=0;
treeFileSave = out2;

doSwapping 			= 1;

originalValueToBeat = valueToBeat;

svl = VERBOSITY_LEVEL;
VERBOSITY_LEVEL = -1;
		    	  
phaseCounter 	= 1;
currentBestTree = treeString;

Tree originalTree = treeString;
treeString = RerootTree (treeString, 0);
Tree givenTree = treeString;

if (MPI_NODE_COUNT>1)
{
	MPINodeState = {MPI_NODE_COUNT-1,1};
	MPINodeTree  = {MPI_NODE_COUNT-1,2};
	MPINodeTree[0]  = "";
}

while (doSwapping)
{
	doSwapping    = 0;
	intBranchCount  = BranchCount ( givenTree );

	fprintf	 (stdout, "\n> num int branches ",intBranchCount,"\n\n"); 

	fprintf	 (stdout, "\n> PHASE ",phaseCounter,"\n\n"); 

	found=0;

	for (ibc = 0; ibc < intBranchCount; ibc=ibc+1)
	{
		fprintf	 (stdout, "\n\n NEW SWAP: \n",givenTree,"\n\n");

		anIntBranch   = BranchName (givenTree,ibc);
		aRerootedTree = RerootTree (givenTree, anIntBranch);
			
		dummy = getTheClustersToSwap (aRerootedTree);
		branchSwap = "(("+clusterOne+","+clusterThree+"),"+clusterTwo+","+clusterFour+")";
		fprintf	 (stdout, "\n FIRST Try: \n",branchSwap,"\n\n");
		dummy = StartSwappingJob(0);														
		if (ibc == intBranchCount)
		{
			break;
		}
		branchSwap = "(("+clusterOne+","+clusterFour+"),"+clusterTwo+","+clusterThree+")";
		fprintf	 (stdout, "\n SECOND Try: \n",branchSwap,"\n\n");
		dummy = StartSwappingJob(0);

	}

	doSwapping=found;
		
	if (doSwapping)
	{
		currentBestTree = RerootTree (BestTreeEver, 0);
		Tree givenTree = currentBestTree;
		fprintf (stdout, "\n\n**** AFTER PHASE ",phaseCounter,"\n\n>Best Tree:\n",BestTreeEver,"\n>Log-likelihood = ",valueToBeat," (Improvement of ", valueToBeat-originalValueToBeat,")\n\n"); 
		phaseCounter = phaseCounter+1;
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





























/*---------------------------------------------------------------*/

function cutTheStringUp		  (delimiterChar)
{
	while (treeStringIndex<treeStringLength)
	{
		aTreeChar = aTreeString[treeStringIndex];
		if ((aTreeChar == delimiterChar)&&(parenthesesDepth==0))
		{
			break;
		}
		else
		{
			if (aTreeChar == ")")
			{
				parenthesesDepth = parenthesesDepth-1;
			}
			else
			{
				if (aTreeChar == "(")
				{
					parenthesesDepth = parenthesesDepth+1;
				}
			}
		}
		treeStringIndex = treeStringIndex+1;
	}
	return 0;
}

/*---------------------------------------------------------------*/

function getTheClustersToSwap (aTreeString)
{
	clusterOne   = "";
	clusterTwo   = "";
	clusterThree = "";
	clusterFour  = "";
	
	treeStringLength = Abs (aTreeString);
	aTreeString		 = aTreeString [2][treeStringLength-3];
	parenthesesDepth = 0;
	treeStringIndex  = 0;
	treeStringLength = treeStringLength - 4;
	startCuttingAt   = 0;
	cutTheStringUp (",");
	clusterOne = aTreeString[startCuttingAt][treeStringIndex-1];
	if (Abs(clusterOne)==0)
	{
		return 1;
	}
	treeStringIndex = treeStringIndex+1;
	startCuttingAt   = treeStringIndex;
	cutTheStringUp (")");
	clusterTwo = aTreeString[startCuttingAt][treeStringIndex-1];
	if (Abs(clusterTwo)==0)
	{
		return 2;
	}
	treeStringIndex = treeStringIndex + 1;
	cutTheStringUp (",");
	treeStringIndex = treeStringIndex+2;
	startCuttingAt   = treeStringIndex;
	cutTheStringUp (",");
	clusterThree = aTreeString[startCuttingAt][treeStringIndex-1];
	if (Abs(clusterThree)==0)
	{
		return 3;
	}
	startCuttingAt = treeStringIndex+1;
	cutTheStringUp (")");
	clusterFour = aTreeString[startCuttingAt][treeStringIndex-1];
	return 4;
}

/*---------------------------------------------------------------*/

function ReceiveJobs (sendOrNot)
{
	if (MPI_NODE_COUNT>1)
	{
		MPIReceive (-1, fromNode, result_String);
		branchSwap_2 = MPINodeTree[fromNode-1][0];
		anIntBranch2 = MPINodeTree[fromNode-1][1];
		if (sendOrNot)
		{
			MPISend (fromNode,theLF);
			MPINodeTree[fromNode-1][0] = branchSwap;
			MPINodeTree[fromNode-1][1] = anIntBranch;
		}
		else
		{
			MPINodeState[fromNode-1]    = 0;
			MPINodeTree[fromNode-1][0]  = "";
			MPINodeTree[fromNode-1][1]  = "";
		}
		
		anIntBranch = anIntBranch2;
		branchSwap	= branchSwap_2;
		ExecuteCommands (result_String);
		
		res = theLF_MLES;
	}
	
	fprintf (stdout, totalRearrangementsTried, "). Swap at:", anIntBranch, ".Log-L = ", Format(res[1][0],10,5)," (", Format(res[1][0]-originalValueToBeat,10,5), ")\n");		
	
	dummy = CheckForImprovement(0);
	
	return fromNode-1;
}

/*---------------------------------------------------------------*/

function CheckForImprovement (dummy)
{
	if (res[1][0]>valueToBeat+0.5*OPTIMIZATION_PRECISION)
	{
		if (globalOption == 2)
		{
			ClearConstraints (testTree);
			ExecuteCommands (setString);
			Optimize (res, theLF);
		}
		totalSupportWeight = totalSupportWeight*Exp(valueToBeat-res[1][0])+1;
				
		valueToBeat     = res[1][0];
		currentBestTree = branchSwap;
		doSwapping 		= 1; 
		savedMLEs		= res;
		
		if (!swapType)
		{
			ibc = intBranchCount;
		}
	}
	else
	{
		totalSupportWeight = totalSupportWeight + Exp(res[1][0]-valueToBeat);
	}
	return 0;
}

/*---------------------------------------------------------------*/

















function StartSwappingJob (dummy)
{
	branchSwap=RerootTree (branchSwap, 0);
	Tree   testTree2 = branchSwap;
	Lea2= TipCount(testTree2);
	Int2= BranchCount(testTree2);
	BranchCount2= Lea2 + Int2;
	fprintf (stdout, "\n Leaves, Internal nodes, total nodes:	\t",Lea2,"\t",Int2,"\t",BranchCount2,"\n");

	for (ibc2 = 0; ibc2 < BranchCount2; ibc2=ibc2+1)											/*iterating over roots*/
	{
		if (ibc2<Int2)
		{
			Tree Tree2 = branchSwap;
			aBranch2   = BranchName (Tree2,ibc2);
			testTreest3 = RerootTree (Tree2, aBranch2);
			Tree testTree=testTreest3;
			fprintf (stdout, "\n Internal node, tree:	",testTree,"\n");
		}
		else
		{
			Tree Tree2 = branchSwap;
			aBranch2   = TipName (Tree2,ibc2-Int2);
			testTreest3 = RerootTree (Tree2, aBranch2);
			Tree testTree=testTreest3;
			fprintf (stdout, "\n External node, tree:	",testTree,"\n");
		}

		LikelihoodFunction theLF = (filteredData, testTree);
		MolecularClock (testTree,t);
		Optimize (res, theLF);
		fprintf (stdout, "\n One root, after optimization:\n",theLF,"\n\n");
		if (res[1][0]>bestabs)
		{
			fprintf (stdout, "\n\n\n FOUND IMPROVEMENT:\n",theLF,"\n\n\n\n\n");
			if (res[1][0]>(bestabs+0.01))
			{
				found=1;
			}
			bestabs=res[1][0];
			valueToBeat=res[1][0];
			fprintf (stdout, "\n\n\n FOUND IMPROVEMENT:\n",theLF,"\n\n\n\n");
			fprintf (treeFileSave,theLF,"\n\n");
			BestTreeEver = testTreest3;
		}
	}
	totalRearrangementsTried = totalRearrangementsTried + 1;
	return 0;
}

