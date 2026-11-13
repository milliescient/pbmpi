#include <cstdlib>

/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include <cassert>
#include "Parallel.h"
#include <string.h>

#include <vector>
#include <sstream>
#include "MultiGeneMixture.h"

// <run name>.<alignment file name without its directory>.treelist
string MultiGeneMixture::GeneTreeFile(int gene)	{
	string base = genename[gene];
	size_t slash = base.find_last_of('/');
	if (slash != string::npos)	{
		base = base.substr(slash + 1);
	}
	return name + "." + base + ".treelist";
}

int MultiGeneMixture::FindGene(const string& s)	{
	for (int gene=0; gene<Ngene; gene++)	{
		string base = genename[gene].substr(genename[gene].find_last_of('/') + 1);
		if ((s == genename[gene]) || (s == base))	{
			return gene;
		}
	}
	if (! myid)	{
		cerr << "warning: " << s << " is not one of the alignments\n";
	}
	return -1;
}

// tree effort per gene, read at every start and restart: <run name>.weights has lines
// "alignment weight" (default 1); a gene in <run name>.frozen has weight 0. A weight w
// gives floor(w) rounds of branch length and SPR moves per cycle, plus one more with
// probability w - floor(w); a gene with no round keeps its tree but still updates its
// mapping and rates.
void MultiGeneMixture::ReadWeights()	{
	genetreeweight = new double[Ngene];
	genefrozen = new int[Ngene];
	for (int gene=0; gene<Ngene; gene++)	{
		genetreeweight[gene] = 1.0;
	}
	ifstream ws((name + ".weights").c_str());
	string s;
	double w;
	while (ws >> s >> w)	{
		int gene = FindGene(s);
		if (gene >= 0)	{
			genetreeweight[gene] = (w > 0) ? w : 0;
		}
	}
	ifstream fs((name + ".frozen").c_str());
	while (fs >> s)	{
		int gene = FindGene(s);
		if (gene >= 0)	{
			genetreeweight[gene] = 0;
		}
	}
	int nfrozen = 0;
	int nweighted = 0;
	for (int gene=0; gene<Ngene; gene++)	{
		genefrozen[gene] = (genetreeweight[gene] == 0);
		nfrozen += genefrozen[gene];
		nweighted += ((genetreeweight[gene] != 0) && (genetreeweight[gene] != 1.0));
	}
	if ((nfrozen || nweighted) && (! myid))	{
		cerr << nfrozen << " genes frozen, " << nweighted << " with a tree weight other than 1, of " << Ngene << '\n';
	}
}

void MultiGeneMixture::MakeFiles()	{
	for (int gene=0; gene<Ngene; gene++)	{
		ofstream os(GeneTreeFile(gene).c_str());
	}
}

void MultiGeneMixture::AllocateAlignments(string datafile, string treefile, int dc)	{

	ifstream is(datafile.c_str());
	is >> Ngene;
	ifstream* tis = 0;
	if (treefile != "None")	{
		tis = new ifstream(treefile.c_str());
		int n;
		(*tis) >> n;
		if (n != Ngene)	{
			cerr << "error: non matching number of alignments and trees\n";
			exit(1);
		}
	}
	genename = new string[Ngene];
	treename = new string[Ngene];
	genesize = new int[Ngene];
	genealloc = new int[Ngene];
	int* geneweight = new int[Ngene];
	// genedata = new SequenceAlignment*[Ngene];
	for (int gene=0; gene<Ngene; gene++)	{
		is >> genename[gene];
		if (tis)	{
			(*tis) >> treename[gene];
		}
		else	{
			treename[gene] = "None";
		}
		SequenceAlignment* data = new FileSequenceAlignment(genename[gene],0,myid);
		if (dc)	{
			data->DeleteConstantSites();
		}
		int nstate = data->GetNstate();
		if (! gene)	{
			Nstate = nstate;
			statespace = data->GetStateSpace();
		}
		else	{
			if (Nstate != nstate)	{
				cerr << "error: all data files do not have the same alphabet\n";
				exit(1);
			}
		}

		genesize[gene] = data->GetNsite();
		// geneweight[gene] = data->GetNsite();
		geneweight[gene] = data->GetNsite() * data->GetNtaxa();
		delete data;
	}
	delete tis;
	tis = 0;

	ReadWeights();
	for (int gene=0; gene<Ngene; gene++)	{
		double w = genetreeweight[gene] > 0.2 ? genetreeweight[gene] : 0.2;
		geneweight[gene] = (int) (geneweight[gene] * w + 0.5);
	}

	// sort alignments by decreasing size
	int permut[Ngene];
	for (int gene=0; gene<Ngene; gene++)	{
		permut[gene] = gene;
	}
	for (int i=0; i<Ngene; i++)	{
		for (int j=Ngene-1; j>i; j--)	{
			if (geneweight[permut[i]] < geneweight[permut[j]])	{
			// if (genesize[permut[i]] < genesize[permut[j]])	{
				int tmp = permut[i];
				permut[i] = permut[j];
				permut[j] = tmp;
			}
		}
	}

	int totsize[nprocs];
	for (int i=0; i<nprocs; i++)	{
		totsize[i] = 0;
	}

	for (int i=0; i<Ngene; i++)	{
		int gene = permut[i];
		int size = geneweight[gene];
		// int size = genesize[gene];

		int min = 0;
		int jmin = 0;
		for (int j=1; j<nprocs; j++)	{
			if ((j==1) || (min > totsize[j]))	{
				min = totsize[j];
				jmin = j;
			}
		}
		genealloc[gene] = jmin;
		totsize[jmin] += size;
	}

	if (totsize[0])	{
		cerr << "error in alloc\n";
		exit(1);
	}
	int total = 0;
	for (int i=1; i<nprocs; i++)	{
		int tot = 0;
		for (int gene=0; gene<Ngene; gene++)	{
			if (genealloc[gene] == i)	{
				tot += geneweight[gene];
				// tot += genesize[gene];
				total++;
			}
		}
		if (tot != totsize[i])	{
			cerr << "error in allocation\n";
			cerr << tot << '\t' << totsize[i] << '\n';
			exit(1);
		}
	}
	if (total != Ngene)	{
		cerr << "error in total allocation\n";
		exit(1);
	}

	genelnL = new double[Ngene];
	tmpgenelnL = new double[Ngene];

	globalnsite = new int[nprocs];
	for (int i=0; i<nprocs; i++)	{
		globalnsite[i] = 0;
	}
	for (int gene=0; gene<Ngene; gene++)	{
		if ((genealloc[gene] < 0) || (genealloc[gene] >= nprocs))	{
			cerr << "alloc : " << genealloc[gene] << '\t' << gene << '\n';
			exit(1);
		}
		globalnsite[0] += genesize[gene];
		globalnsite[genealloc[gene]] += genesize[gene];
	}
	if (! myid)	{
		cerr << '\n';
		for (int i=1; i<nprocs; i++)	{
			cerr << i << '\t' << globalnsite[i] << '\n';
		}
		cerr << '\n';
		cerr << "total: " << GetGlobalNsite() << '\n';
	}
	
	// check total size
	if (! myid)	{
		int tot = 0;
		for (int i=1; i<nprocs; i++)	{
			tot += globalnsite[i];
		}
		if (tot != globalnsite[0])	{
			cerr << "error in total size during gene allocation\n";
			exit(1);
		}
		int tot2 = 0;
		for (int gene=0; gene<Ngene; gene++)	{
			tot2 += genesize[gene];
		}
		if (tot2 != tot)	{
			cerr << "error during alloc: total size does not match\n";
			exit(1);
		}
	}

	genealpha = new double[Ngene];
	genelength = new double[Ngene];
	delete[] geneweight;
}

double MultiGeneMixture::GetLogLikelihood()	{
	GlobalCollectGeneLikelihoods();
	lnL = 0;
	for (int gene=0; gene<Ngene; gene++)	{
		lnL += genelnL[gene];
	}
	return lnL;
}

void MultiGeneMixture::GlobalCollectGeneLikelihoods()	{
	// send signal
	assert(myid == 0);
	MESSAGE signal = LIKELIHOOD;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

	MPI_Status stat;
	for (int j=1; j<nprocs; j++)	{
		MPI_Recv(tmpgenelnL,Ngene,MPI_DOUBLE,j,TAG1,MPI_COMM_WORLD,&stat);
		for (int gene=0; gene<Ngene; gene++)	{
			if (genealloc[gene] == j)	{
				genelnL[gene] = tmpgenelnL[gene];
			}
		}
	}
}

void MultiGeneMixture::SlaveSendGeneLikelihoods()	{
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			genelnL[gene] = process[gene]->GetLogLikelihood();
		}
	}
	MPI_Send(genelnL,Ngene,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
}

double MultiGeneMixture::GetMeanLength()	{
	GlobalCollectGeneLengths();
	double total = 0;
	for (int gene=0; gene<Ngene; gene++)	{
		total += genelength[gene];
	}
	return total / Ngene;
}

void MultiGeneMixture::GlobalCollectGeneLengths()	{
	// send signal
	assert(myid == 0);
	MESSAGE signal = LENGTH;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

	MPI_Status stat;
	for (int j=1; j<nprocs; j++)	{
		MPI_Recv(tmpgenelnL,Ngene,MPI_DOUBLE,j,TAG1,MPI_COMM_WORLD,&stat);
		for (int gene=0; gene<Ngene; gene++)	{
			if (genealloc[gene] == j)	{
				genelength[gene] = tmpgenelnL[gene];
			}
		}
	}
}

void MultiGeneMixture::SlaveSendGeneLengths()	{
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			genelength[gene] = process[gene]->GetTotalLength();
		}
	}
	MPI_Send(genelength,Ngene,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
}

double MultiGeneMixture::GetMeanAlpha()	{
	GlobalCollectGeneAlphas();
	double total = 0;
	for (int gene=0; gene<Ngene; gene++)	{
		total += genealpha[gene];
	}
	return total / Ngene;
}

void MultiGeneMixture::GlobalCollectGeneAlphas()	{
	// send signal
	assert(myid == 0);
	MESSAGE signal = ALPHA;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

	MPI_Status stat;
	for (int j=1; j<nprocs; j++)	{
		MPI_Recv(tmpgenelnL,Ngene,MPI_DOUBLE,j,TAG1,MPI_COMM_WORLD,&stat);
		for (int gene=0; gene<Ngene; gene++)	{
			if (genealloc[gene] == j)	{
				genealpha[gene] = tmpgenelnL[gene];
			}
		}
	}
}

void MultiGeneMixture::SlaveSendGeneAlphas()	{
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			genealpha[gene] = process[gene]->GetAlpha();
		}
	}
	MPI_Send(genealpha,Ngene,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
}

void MultiGeneMixture::GlobalGeneMove()	{
	assert(myid == 0);
	MESSAGE signal = GENE_MOVE;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
}

void MultiGeneMixture::GlobalSample()	{
	assert(myid == 0);
	MESSAGE signal = SAMPLE;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
}

void MultiGeneMixture::GlobalUnfold()	{
	assert(myid == 0);
	MESSAGE signal = UNFOLD;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
}

void MultiGeneMixture::GlobalCollapse()	{
	assert(myid == 0);
	MESSAGE signal = COLLAPSE;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
}

void MultiGeneMixture::SaveTrees()	{
	assert(myid == 0);
	MESSAGE signal = SAVETREES;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
}

void MultiGeneMixture::GlobalFromStream(istream& is)	{
	assert(myid == 0);
	MESSAGE signal = FROMSTREAM;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

	string tag;
	is >> tag;
	if (tag == "MULTIGENE_STATE")	{
		vector<string> state(Ngene);
		for (int k=0; k<Ngene; k++)	{
			int gene;
			unsigned int glen;
			is >> gene >> glen;
			is.get();
			string g(glen,' ');
			is.read(&g[0],glen);
			state[gene] = g;
		}
		for (int j=1; j<nprocs; j++)	{
			ostringstream os;
			os << "MGREC\n";
			for (int gene=0; gene<Ngene; gene++)	{
				if (genealloc[gene] == j)	{
					os << gene << '\t' << state[gene].size() << '\n' << state[gene];
				}
			}
			string s = os.str();
			unsigned int len = s.length();
			MPI_Send(&len,1,MPI_INT,j,TAG1,MPI_COMM_WORLD);
			MPI_Send(&s[0],len,MPI_UNSIGNED_CHAR,j,TAG1,MPI_COMM_WORLD);
		}
	}
	else	{
		// earlier format, one block per worker: needs the same workers and allocation
		for (int j=1; j<nprocs; j++)	{
			unsigned int len = 0;
			if (j == 1)	{
				len = atoi(tag.c_str());
			}
			else	{
				is >> len;
			}
			string s(len,' ');
			for (unsigned int i=0; i<len; i++)	{
				is.get(s[i]);
			}
			MPI_Send(&len,1,MPI_INT,j,TAG1,MPI_COMM_WORLD);
			MPI_Send(&s[0],len,MPI_UNSIGNED_CHAR,j,TAG1,MPI_COMM_WORLD);
		}
	}
}

void MultiGeneMixture::SlaveFromStream()	{

	MPI_Status stat;
	unsigned int len = 0;
	MPI_Recv(&len, 1, MPI_INT,0,TAG1,MPI_COMM_WORLD,&stat);
	string s(len,' ');
	MPI_Recv(&s[0],len,MPI_UNSIGNED_CHAR,0,TAG1,MPI_COMM_WORLD,&stat);
	if (s.compare(0,6,"MGREC\n") == 0)	{
		istringstream bs(s.substr(6));
		int gene;
		unsigned int glen;
		while (bs >> gene >> glen)	{
			bs.get();
			string g(glen,' ');
			bs.read(&g[0],glen);
			istringstream gs(g);
			process[gene]->FromStream(gs);
		}
	}
	else	{
		istringstream is(s);
		for (int gene=0; gene<Ngene; gene++)	{
			if (genealloc[gene] == myid)	{
				process[gene]->FromStream(is);
			}
		}
	}
}

// gene states are saved one record per gene, by gene index, so that a restart can give
// genes to different workers: a line MULTIGENE_STATE, then for each gene its index, the
// length of its state, a newline and the state itself
void MultiGeneMixture::GlobalToStream(ostream& os)	{

	assert(myid == 0);
	MESSAGE signal = TOSTREAM;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

	vector<string> state(Ngene);
	MPI_Status stat;
	for (int j=1; j<nprocs; j++)	{
		unsigned int len = 0;
		MPI_Recv(&len, 1, MPI_INT,j,TAG1,MPI_COMM_WORLD,&stat);
		string blob(len,' ');
		MPI_Recv(&blob[0],len,MPI_UNSIGNED_CHAR,j,TAG1,MPI_COMM_WORLD,&stat);
		istringstream bs(blob);
		int gene;
		unsigned int glen;
		while (bs >> gene >> glen)	{
			bs.get();
			string g(glen,' ');
			bs.read(&g[0],glen);
			state[gene] = g;
		}
	}
	os << "MULTIGENE_STATE\n";
	for (int gene=0; gene<Ngene; gene++)	{
		os << gene << '\t' << state[gene].size() << '\n' << state[gene];
	}
}

void MultiGeneMixture::SlaveToStream()	{
	ostringstream os;
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			ostringstream gs;
			process[gene]->ToStream(gs);
			string g = gs.str();
			os << gene << '\t' << g.size() << '\n' << g;
		}
	}
	string s = os.str();
	unsigned int len = s.length();
	MPI_Send(&len,1,MPI_INT,0,TAG1,MPI_COMM_WORLD);
	MPI_Send(&s[0],len,MPI_UNSIGNED_CHAR,0,TAG1,MPI_COMM_WORLD);
}

void MultiGeneMixture::SlaveGeneMove()	{
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			double w = genetreeweight[gene];
			int n = (int) w;
			if (rnd::GetRandom().Uniform() < w - n)	{
				n++;
			}
			if (n)	{
				process[gene]->SetTreeReps(n);
				process[gene]->Move();
			}
			else	{
				process[gene]->FrozenMove();
			}
		}
	}
}

void MultiGeneMixture::SlaveSample()	{
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->Sample();
		}
	}
}

void MultiGeneMixture::SlaveUnfold()	{
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->Unfold();
		}
	}
}

void MultiGeneMixture::SlaveCollapse()	{
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->Collapse();
		}
	}
}

void MultiGeneMixture::SlaveSaveTrees()	{
	for (int gene=0; gene<Ngene; gene++)	{
		if ((genealloc[gene] == myid) && (! genefrozen[gene]))	{
			process[gene]->SaveTree(GeneTreeFile(gene));
		}
	}
}

void MultiGeneMixture::WaitLoop()	{
	MESSAGE signal;
	do {
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
		if (getenv("MPB_TRACE_SIGNALS")) cerr << "slave " << myid << " signal " << signal << '\n';
		if (signal == KILL) break;
		SlaveExecute(signal);
	} while(true);
}

void MultiGeneMixture::SlaveExecute(MESSAGE signal)	{

	switch(signal) {
	case FROMSTREAM:
		SlaveFromStream();
		break;
	case TOSTREAM:
		SlaveToStream();
		break;
	case SAVETREES:
		SlaveSaveTrees();
		break;
	case LIKELIHOOD:
		SlaveSendGeneLikelihoods();
		break;
	case LENGTH:
		SlaveSendGeneLengths();
		break;
	case ALPHA:
		SlaveSendGeneAlphas();
		break;
	case SAMPLE:
		SlaveSample();
		break;
	case UNFOLD:
		SlaveUnfold();
		break;
	case COLLAPSE:
		SlaveCollapse();
		break;
	case UPDATE_SPROFILE:
		SlaveUpdateSiteProfileSuffStat();
		break;
	case PARAMETER_DIFFUSION:
		SlaveUpdateParameters();
		break;
	case GENE_MOVE:
		SlaveGeneMove();
		break;
	case MIX_MOVE:
		SlaveMixMove();
		break;
	
	default:
		cerr << "slave " << myid << ": unhandled signal " << signal << '\n';
		exit(1);
	}
}

