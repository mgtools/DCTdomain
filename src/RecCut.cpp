#include <bits/stdc++.h>
#include <utility>
#include <fstream>
#include <sstream>
#include <vector>
using namespace std;
//Haixu Tang, hatang@indiana.edu

//#define cuts_threshold 0.04
#define Min_Terminal 10

double cuts_threshold=0.08;
double cuts_2_threshold=0.07;
int Min_Size=22;

void SplitDomain(int cuts, std::vector<std::pair<int, int>>& domains, std::vector<int>& cutsite, std::vector<std::vector<std::pair<int, int>>>& domains1,
	       	std::vector<std::vector<std::pair<int, int>>>& domains2, std::vector<int>& cutsite1, std::vector<int>& cutsite2, int length)
{
	std::vector<std::pair<int, int>> domain1, domain2;
	int i = 0;
	int j = 0;
	int c1;

	for(c1 = 0; c1 < cutsite.size(); c1 ++)	{
		if(cutsite[c1] < cuts)	{
			domain1.push_back(domains[c1]);
			cutsite1.push_back(cutsite[c1]);
		} else {
			break;
		}
	}

	//split the segment c1
	int len1, len2;
	if(c1 == 0)	{
		len1 = cuts;
	} else	{
		len1 = cuts - cutsite[c1 - 1];
	}
	domain1.push_back(std::make_pair(domains[c1].first, domains[c1].first+len1-1));

	// second half of the segment c1 goes to domain 2
	if(cutsite.size() == 0 || c1 == cutsite.size())	{
		len2 = length - cuts;
	} else	{
		len2 = cutsite[c1] - cuts;
	}
	if(len2 > 0)	{
		domain2.push_back(std::make_pair(domains[c1].second-len2+1, domains[c1].second));
		if(c1 < cutsite.size())	{
			cutsite2.push_back(len2);
		}
	}
	for(c1 ++; c1 < domains.size(); c1 ++)	{
		domain2.push_back(domains[c1]);
		if(c1 < cutsite.size() && len1 > 0)	{
			cutsite2.push_back(len1);
		}
	}

	domains1.push_back(domain1);
	domains2.push_back(domain2);
}

void SplitDomain_2cuts(int cuts1, int cuts2, std::vector<std::pair<int, int>>& domains, std::vector<int>& cutsite,
	       	std::vector<std::vector<std::pair<int, int>>>& domains1, std::vector<std::vector<std::pair<int, int>>>& domains2,
	       	std::vector<int>& cutsite1, std::vector<int>& cutsite2, int length)
{
	std::vector<std::pair<int, int>> domain1, domain2;
	int i = 0;
	int j = 0;

	int len = length - cuts2;
	int c1, c2;
	c1 = -1;
	c2 = 0;
	for(i = cutsite.size() - 1; i >= 0; i --)	{
		//printf("cutsite %d %d\n", i, cutsite[i]);
		if(cutsite[i] < cuts2)	{
			c2 = i + 1;
		}
	       	if(cutsite[i] <= cuts1)	{
			c1 = i;
			break;
		}
	}
	// cuts1 before c1
	
	//split the segment (cuts1, cuts2)
	int len1, len2;
	if(c2 == 0)	{
		len2 = cuts2;
	} else	{
		len2 = cuts2 - cutsite[c2 - 1];
	}
	// create the segments in the first domain
	domain1.push_back(std::make_pair(domains[c2].first + len2, domains[c2].second));
	for(i = c2 + 1; i <= cutsite.size(); i ++)	{	// # domains is one above the # cutsitef
		domain1.push_back(domains[i]);
	}
	// copy cuts in (cuts2, V) into cutsite1
	for(i = c2; i < cutsite.size(); i ++)	{
		cutsite1.push_back(cutsite[i] - cuts2);
	}
	// create a new cut between the last segment and the first segment
	cutsite1.push_back(len);
		
	//copy cuts in (0, cuts1) into cutsite 1
	for(i = 0; i <= c1; i ++)	{
		domain1.push_back(domains[i]);
		cutsite1.push_back(cutsite[i] + len);
	}
	if(c1 >= 0)	{
		len1 = cuts1 - cutsite[c1];
	} else	{
		len1 = cuts1;
	}
	c1 ++; // cuts1 within the segment of c1
	if(len1 > 0) {
		domain1.push_back(std::make_pair(domains[c1].first, domains[c1].first + len1 - 1));
	}

	// create the segments in the 2nd domain
	if(c2 > c1 && cutsite[c1 - 1] + len1 < cutsite[c1])	{
		if(c1 >= 0)	{
			domain2.push_back(std::make_pair(domains[c1].first + len1, domains[c1].second));
		}
	}
	for(i = c1 + 1; i < c2; i ++)	{
		domain2.push_back(domains[i]);
	}
	if(c2 > c1 && len2 > 0)	{
		domain2.push_back(std::make_pair(domains[c2].first, domains[c2].first + len2 - 1));
	} else if(c1 == c2)	{ 	// both cuts are in the same segment
		domain2.push_back(std::make_pair(domains[c2].first + len1, domains[c2].first + len2 - 1));
	}
	// copy cuts in (cuts1, cuts2) into cutsite 2
	for(i = c1; i < c2; i ++)	{
		if(cutsite[i] - cuts1 > 0)	{
			cutsite2.push_back(cutsite[i] - cuts1);
		}
	}

	domains1.push_back(domain1);
	domains2.push_back(domain2);
}

void recursiveMaxCut(std::vector<std::vector<int>>& graph, std::vector<int>& cutsite, std::vector<std::vector<std::pair<int, int>>>& domains)
{
    int V = graph.size();
    if(V < Min_Size)	{
	    cout << "Protein has length of 0" << endl;
	    exit(-1);
    }
    int num_domains = domains.size();
    if(num_domains != 1)	{
	   cout << "Error: Number of domains is not 1: " << num_domains << endl;
	   exit(-1);
    }

    double mincut = 2.0;
    int cutv = 0;
    int pre = 0;
    int post = 0;
    int cuts;
    double avecut;
    int	sum = 0;
    for (int i = 0; i < V - 1; i++) {
        for (int j = i + 1; j < V; j++) {
		sum += graph[i][j];
        }
    }

    std::vector<int> N1(V);
    std::vector<int> N2(V);
    N1[0] = 0;
    N2[0] = sum;
    for(int j = 1; j < V; j ++)	{
	    N2[j] -= graph[0][j];
	    cutv += graph[0][j];
    }
    for(int i = 1; i < V - 2; i ++)	{	// move vertex i from the second partition to the 1st partiion
	pre = 0; 
        for(int j = 0; j < i; j ++)	{
	     pre += graph[i][j];
	}
	post = 0; 
        for(int j = i + 1; j < V; j ++)	{
	     post += graph[i][j];
	}
	cutv = cutv + post - pre;
	N1[i] = N1[i - 1] + pre; 
	N2[i] = N2[i - 1] - post;
        avecut = ((double) cutv) * sum / N1[i] / N2[i];
	if(avecut < mincut)	{
		mincut = avecut;
		cuts = i + 1;
	}
    }
    for(int i = V - 2; i < V; i ++)	{	// move vertex i from the second partition to the 1st partiion
	pre = 0; 
        for(int j = 0; j < i; j ++)	{
	     pre += graph[i][j];
	}
	post = 0; 
        for(int j = i + 1; j < V; j ++)	{
	     post += graph[i][j];
	}
	N1[i] = N1[i - 1] + pre; 
	N2[i] = N2[i - 1] - post;
    }

    // Compute scores for double-cuts
    std::vector<std::vector<int>> T(V, std::vector<int>(V, 0));
    std::vector<std::vector<int>> TS(V, std::vector<int>(V, 0));
    std::vector<std::vector<int>> S(V, std::vector<int>(V, 0));
    for(int i = 0; i < V; i ++)	{
	T[i][0] = 0;
    	for(int j = 1; j < V; j ++)	{
		if(j != i)
			T[i][j] = T[i][j - 1] + graph[i][j];
		else
			T[i][j] = T[i][j - 1];
	}
    	for(int j = i + 1; j < V; j ++)	{
		if(i == 0)
			TS[i][j] = T[i][j];
		else
			TS[i][j] = T[i][j] - T[i][i - 1];
	}
    }
    // Compute the sum of contact scores between i and j
    for(int j = 0; j < V; j ++)	{
	S[j][j] = 0;
    	for(int i = j - 1; i >= 0; i --)	{
		S[i][j] = S[i + 1][j] + TS[i][j];
	}
    }

    // Compute double-cut at i and j
    int NS1, NS2;
    int cuts1, cuts2;
    double mincut2;
    mincut2 = 2.0;
    for(int i = Min_Terminal; i < V - Min_Terminal; i ++)	{
    	for(int j = i + Min_Size - 1; j < V - Min_Terminal; j ++)	{
		//NS1 = N1[i] + N2[j + 1];
		NS2 = S[i][j];
		cutv = N1[j] + N2[i] - N1[i - 1] - N2[j + 1] - NS2 * 2;
		NS1 = sum - cutv - NS2;
		if(NS1 > 0 && NS2 > 0)	{
		        avecut = ((double) cutv) * sum / NS1 / NS2;
			if(avecut < mincut2)	{
				mincut2 = avecut;
				cuts1 = i;
				cuts2 = j;
			}
		}
	}
    }

    if(mincut - cuts_threshold <= mincut2 - cuts_2_threshold)	{
       if(mincut > cuts_threshold || cuts < Min_Size || V - cuts < Min_Size)	{
	    	return;
       } else	{
    // Create the 2D vectors for the subgraphs
		std::vector<std::vector<int>> firstgraph(cuts, std::vector<int>(cuts));
		std::vector<std::vector<int>> secondgraph(V-cuts, std::vector<int>(V-cuts));

    // Copy the first cuts x cuts subgraph
                for (int i = 0; i < cuts; i++) {
        		for (int j = 0; j < cuts; j++) {
            			firstgraph[i][j] = graph[i][j];
        		}
    		}
    // Copy the (V - cuts) x (V - cuts) subgraph starting from cuts
    		for (int i = cuts; i < V; i++) {
        		for (int j = cuts; j < V; j++) {
            			secondgraph[i - cuts][j - cuts] = graph[i][j];
        		}
    		}

    		std::vector<int> cutsite1, cutsite2;
		std::vector<std::vector<std::pair<int, int>>> domains1;
		std::vector<std::vector<std::pair<int, int>>> domains2;

		SplitDomain(cuts, domains[0], cutsite, domains1, domains2, cutsite1, cutsite2, V);
	    	recursiveMaxCut(firstgraph, cutsite1, domains1);
	    	recursiveMaxCut(secondgraph, cutsite2, domains2);
		//cout << "num_domains2 " << domains2.size() << endl;
		// Merges domains1 and domains2
		domains = domains1;
		domains.insert(domains.end(), domains2.begin(), domains2.end());
       }
       return;
    } else {
       int length = cuts2 - cuts1;
       if(mincut2 > cuts_2_threshold || length < Min_Size || V - length < Min_Size)	{
	    	return;
       } else	{
    // Create the 2D vectors for the subgraphs
		std::vector<std::vector<int>> firstgraph(V - length, std::vector<int>(V - length));
		std::vector<std::vector<int>> secondgraph(length, std::vector<int>(length));

    // Copy the first subgraph (cuts2, V) and (0, cuts1)	first graph: (cuts2, V) --> (0, cuts1 - 1)
    		int last_size = V - cuts2;	// length of the segment of (cuts2, V)
    		for (int i = 0; i < cuts1; i++) {
        		for (int j = i + 1; j < cuts1; j++) {
            			firstgraph[j + last_size][i + last_size] = firstgraph[i + last_size][j + last_size] = graph[i][j];
				//printf("last_size %d coord %d %d %d %d\n", last_size, i, j,  i + last_size, j + last_size);
				//getchar();
        		}
    		}
    		for (int i = cuts2; i < V; i++) {
        		for (int j = i + 1; j < V; j++) {
            			firstgraph[j - cuts2][i - cuts2] = firstgraph[i - cuts2][j - cuts2] = graph[i][j];
        		}
    		}
    		for (int i = 0; i < cuts1; i++) {
        		for (int j = cuts2; j < V; j++) {
            			firstgraph[j - cuts2][i + last_size] = firstgraph[i + last_size][j - cuts2] = graph[i][j];
        		}
    		}
    // Copy the second subgraph between cuts1 and cuts2
                for (int i = cuts1; i < cuts2; i++) {
        		for (int j = i + 1; j < cuts2; j++) {
            			secondgraph[j - cuts1][i - cuts1] = secondgraph[i - cuts1][j - cuts1] = graph[i][j];
        		}
    		}

    		std::vector<int> cutsite1, cutsite2;
		std::vector<std::vector<std::pair<int, int>>> domains1;
		std::vector<std::vector<std::pair<int, int>>> domains2;

		SplitDomain_2cuts(cuts1, cuts2, domains[0], cutsite, domains1, domains2, cutsite1, cutsite2, V);
		//output two domains
		int V1 = firstgraph.size();
		//printf("Graph 1 size: %d\n", V1);
	    	recursiveMaxCut(firstgraph, cutsite1, domains1);
		int V2 = secondgraph.size();
		//printf("Graph 2 size: %d\n", V2);
	    	recursiveMaxCut(secondgraph, cutsite2, domains2);
		// Merges domains1 and domains2
		domains = domains1;
		domains.insert(domains.end(), domains2.begin(), domains2.end());
       }
       return;
    }
}


std::vector<std::vector<int>> readGraph(const std::string& filename) {
    std::ifstream file(filename);
    
    // Read the first line and parse the number of vertices
    std::string line;
    std::getline(file, line);
    std::stringstream ss(line);
    std::string item, item1;
    int num_vertices=0;
    ss >> item >> item1 >> num_vertices; // Read items separated by spaces
    if(num_vertices == 0)	{
	    cout << "No protein input" << endl;
	    exit(-1);
    }

    // Create a 2D vector with all elements initialized to 0 (not contact)
    std::vector<std::vector<int>> graph(num_vertices, std::vector<int>(num_vertices, 0));
    int V = graph.size();

    int count = 0;
    // Parse the file line by line
    while (std::getline(file, line)) {
        // If this line starts with "CON", parse the graph data
        if (line.substr(0,3) == "CON") {
            std::stringstream ss(line.substr(3));
            std::string item;
            while (std::getline(ss, item, ',')) {
                int i, j;
		double value;
                std::stringstream(item) >> i >> j >> value;
	        graph[i][j] = graph[j][i] = (int) (value * 100 + 0.5);
		count ++;
            }
        }
    }
    for(int i = 0; i < V; i ++)	{
	//for(int j = 1; j <= 4 && i+j < V; j ++)	{
	for(int j = 1; j <= 3 && i+j < V; j ++)	{
	    graph[i][i + j] = graph[i + j][i] = 100;	//adjacent residues must be in contact
	}
    }

    return graph;
}


int main(int argc, char *argv[])
{
    int i;
    int filename_given = 0;
    int showname_given = 0;
    std::string filename, showname;
    for(i = 0; i < argc; i ++) {
        if((!strcmp(argv[i], "--input")) && (argc > i + 1)) {
	   filename = argv[i + 1];
	   filename_given = 1;
	}
	else if((!strcmp(argv[i], "--name")) && (argc > i + 1)) {
	   showname = argv[i + 1];
	   showname_given = 1;
	}
	else if((!strcmp(argv[i], "--cutoff")) && (argc > i + 2)) {
	   cuts_threshold = atof(argv[i + 1]);
	   cuts_2_threshold = atof(argv[i + 2]);
        }
    }
    if (!(filename_given && showname_given)) {
        std::cerr << "Usage: <--input input-CE-file> <--name name-to-show> [--cutoff continuous_domains discontinuous_domains]\n";
        return 1;
    }

    std::vector<std::vector<int>> graph = readGraph(filename);
    int V = graph.size();
    //std::cout << "# vertices: " << V << endl;

    // Initially, all vertices are included.
    std::vector<int> vertices(V);
    std::iota(vertices.begin(), vertices.end(), 0);

    // Initialize the partition indicators. Start with all 0s.
    std::vector<int> partitions(V, 0);

    // Start the recursive process with partition 1.
    //int total_partition;
    //total_partition = recursiveMinCut(graph, vertices, partitions, 1);
    std::vector<int> cutsite;
    std::vector<std::vector<std::pair<int, int>>> domains;
    vector<std::pair<int, int>> domain;

    domain.push_back(std::make_pair(0, V-1));
    domains.push_back(domain);
    if(graph.size() >= Min_Size)    
    	recursiveMaxCut(graph, cutsite, domains);

    //output the domains: 1-70,180-261,426-511; 262-425; 71-179;
    std::cout<<showname<<" "<<domains.size()<<" ";
    pair<int, int> p;
    for (i = 0; i <  domains.size(); i ++) {
	domain = domains[i];
	for(int j = 0; j < domain.size(); j ++)	{
		p = domain[j];
		if(j > 0)	std::cout << ",";
	    	std::cout << p.first + 1<< "-" << p.second + 1;
	}
	std::cout << ";";
    }
    std::cout << endl;

    return 0;
}
