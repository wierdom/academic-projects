#include <cassert>

#include <cstdlib>
using std::exit;

#include <cstddef>
using std::size_t;
using std::ptrdiff_t;

#include <algorithm>
using std::min;
using std::max;
using std::fill;
using std::shuffle;
using std::random_shuffle;

#include <utility>
using std::pair;
using std::make_pair;

#include <iostream>
using std::endl;
using std::cout;
using std::cerr;

#include <fstream>
using std::ifstream;
using std::ofstream;
	
#include <iomanip>
using std::setw;

#include <set>
using std::set;

#include <queue>
using std::queue;

#include <random>
using std::mt19937;
using std::uniform_real_distribution;

#include <chrono>
using std::chrono::high_resolution_clock;

#include <functional>
using std::bind;

// signal handling
#include <csignal>
#include <unistd.h>

// chain
//const size_t Lx = 64;   // Need an even number of physical sites
//const size_t Ns = Lx;   // Even number of physical sites

// square
const size_t Lx = 8; // Need an even number of physical sites along each direction
const size_t Ly = 8; // Need an even number of physical sites along each direction
const size_t Ns = Lx*Ly;   // Doubly-even number of physical sites

const size_t twiceS = 4;      // (Integer) number of bonds emerging from each site
const size_t D = Ns*twiceS;   // Number of effective S=1/2 degrees of freedom

//   0     1     2     3     ...   Ns
//   ----  ----  ----  ----        ----
//   n=0   2S    4S    6S    ...   2S*(Ns-1)
//   1     2S+1  4S+1  6S+1        .
//   2     2S+2  4S+2  .           .
//   .     .     .     .
//   .     .     .
//   2S-1  4S-1  6S-1              2S*Ns-1

static size_t left[D];
static size_t right[D];

static int label[D];
static int ordering[D];
static int pseudospin[D];
static int sublat[D]; // specify sublattice (0 or 1)

static int num_loops;
static queue<int> unused_labels;

static size_t site_index[Ns];
const size_t num_flavour_index_pairs = twiceS*(twiceS-1)/2;
static pair<size_t,size_t> flavour_index_pairs[num_flavour_index_pairs];

static double spin_corr[Ns];
static double spin_corr_collect[Ns];

static mt19937 rng; // mersenne twistor random number generator
static auto rnd = bind(uniform_real_distribution<double>(0,1),rng); // random number in [0,1)

void init_random()
{
	ifstream frand("mtrand_internal_state.dat");
	if (frand.is_open())
		frand >> rng;
	else
		rng.seed(mt19937::result_type(high_resolution_clock::now().time_since_epoch().count()));
	frand.close();
}

void shut_down(void)
{
	// save internal state of the random number generator
	ofstream frand("mtrand_internal_state.dat");
	if (frand.is_open())
		frand << rng;
	frand.close();
	exit(1);
}

void interrupt_handler(int)
{
	shut_down();
}

void init_interrupt_handler(void)
{
	// prepare to catch ctrl-c events
	struct sigaction SIGINT_handler;
	SIGINT_handler.sa_handler = interrupt_handler;
	sigemptyset(&SIGINT_handler.sa_mask);
	SIGINT_handler.sa_flags = 0;
	sigaction(SIGINT, &SIGINT_handler, NULL);
}

void init_config_triplet(size_t bond[D])
{
	if (twiceS%2 == 0) // only valid for even number of spins-1/2 per site
	{
		assert(D%(twiceS)==0);
		// generic onsite triplet state for integer S
		//  (0S)~~(0S+1)~ ... ~(2S-1)
		//  (2S)~~(2S+1)~ ... ~(4S-1)
		//   .       .            .
		//   .       .            .
		//   .       .            .
		//  (D-2S)~~(D-2S+1)~ ... ~(D-1)
		for (size_t i = 0; i < Ns; ++i)
		{
			for (size_t f = 0; f < twiceS/2; ++f)
			{
				size_t a=i*twiceS+2*f;
				size_t b=i*twiceS+2*f+1;
				sublat[a] = i%2;
				sublat[b] = i%2;
				bond[a] = b;
				bond[b] = a;
				//cout << "a=" << a << " b=" << b << endl;
			}
		}
	}
	else 
	{
		cerr << "2S = " << twiceS 
			 << " is not a valid case in function init_config_triplet_chain()." << endl;
		exit(1);
	}
}

void init_config_aklt_square(size_t bond[D])
{
        assert(Lx > 3); // minimum number of sites along x direction
        assert(Lx%2 == 0); // even number of sites along x direction
        assert(Ns%4 == 0);
	if (twiceS%4 == 0) // 2S a multiple of z=4
	{
		// specialized to the square lattice geometry
		//   !   !   !   !
		// ~~00~~01~~02~~03~~
		//   !   !   !   !
		// ~~07~~06~~05~~04~~
		//   !   !   !   !
		// ~~08~~09~~10~~11~~
		//   !   !   !   !
		// ~~15~~14~~13~~12~~
		//   !   !   !   !
		for (size_t n = 0; n < Ns; ++n)
		{
			size_t h=(n+1);
			if (h%Lx == 0) h=h-Lx;
			size_t v=(n+2*(Lx-n%Lx)-1)%Ns;
			//cout << "n=" << n << " h=" << h << " v=" << v << endl;
			for (size_t f = 0; f < twiceS/2; ++f)
			{
				size_t a=2*f;
				size_t b=2*f+1;
				a+=n*twiceS;
				if (f%2 == 0)
				{
					b+=h*twiceS;
					bond[a] = b; // horizontal
					bond[b] = a; // horizontal
				}
				else
				{
					b+=v*twiceS;
					bond[a] = b; // vertical
					bond[b] = a; // vertical
				}
			}
			for (size_t f = 0; f < twiceS; ++f)
			{
				size_t a=f+n*twiceS;
				sublat[a] = n%2;
			}
		}
	}
	else 
	{
		cerr << "2S = " << twiceS 
			 << " is not a valid case in function init_config_aklt_square()." << endl;
		exit(1);
	}
}

void verify_config(size_t bond[D])
{
	set<size_t> S;
	for (size_t n = 0; n < D; ++n)
	{
		assert(bond[n] != n);
		assert(bond[n] >= 0);
		assert(bond[n] <= D);
		assert(bond[bond[n]] == n);
		S.insert(bond[n]);
	}
	auto step = S.cbegin();
	for (size_t n = 0; n < D; ++n)
	{
		assert(step != S.cend());
		assert(*step == n);
		++step;
	}
	assert(step == S.cend());
}

void write_labels_and_ordering(int o, size_t n, int l)
{
	//int o = 0; // ordering around each loop
	size_t walk = n;
	do
	{
		label[walk] = l;
		ordering[walk] = o++;
		walk = left[walk];
		label[walk] = l;
		ordering[walk] = o++;
		walk = right[walk];
	} while (walk != n);
}

void write_ordering(size_t n)
{
	int o = 0; // ordering around each loop
	size_t walk = n;
	do
	{
		ordering[walk] = o++;
		walk = left[walk];
		ordering[walk] = o++;
		walk = right[walk];
	} while (walk != n);
}

void init_labels(void)
{
	fill(ordering,ordering+D,-1); // -1 as a placeholder for `as yet unassigned`
	fill(label,label+D,-1);

	int l = 0; // current label	
	for (size_t n = 0; n < D; ++n)
	{
		if (label[n] == -1) // if unassigned
		{
			assert(ordering[n] == -1); // also unassigned
			write_labels_and_ordering(0,n,l);
			assert(ordering[n] == 0 || ordering[n] == 1);
			assert(label[n] == l);
			++l;
		}
	}
	num_loops = l; // set the total number of loops
	assert(D%2 == 0);
	for ( ; l < int(D)/2; ++l)
		unused_labels.push(l);
}

void init_labels_fast(void)
{
	fill(ordering,ordering+D,-1); // -1 as a placeholder for `as yet unassigned`
	fill(label,label+D,-1);

	int l = 0; // current label	
	for (size_t n = 0; n < D; ++n)
	{
		if (label[n] == -1) // if unassigned
		{
			assert(ordering[n] == -1); // also unassigned
			write_labels_and_ordering(0,n,l);
			assert(ordering[n] == 0 || ordering[n] == 1);
			assert(label[n] == l);
			++l;
		}
	}
	num_loops = l; // set the total number of loops
}

void write_labels_and_pseudospin(int o, size_t n, int l)
{
	size_t walk = n;
	do
	{
		label[walk] = l;
		pseudospin[walk] = o+int(sublat[walk]);
		walk = left[walk];
		label[walk] = l;
		pseudospin[walk] = o+int(sublat[walk]);
		walk = right[walk];
	} while (walk != n);
}

void init_pseudospin(void)
{
	fill(pseudospin,pseudospin+D,-1); // -1 as a placeholder for `as yet unassigned`
	fill(label,label+D,-1);

	int l = 0; // current label	
	for (size_t n = 0; n < D; ++n)
	{
		if (label[n] == -1) // if unassigned
		{
			assert(pseudospin[n] == -1); // also unassigned
			int r=0;
			if (rnd() < 0.5) ++r;
			write_labels_and_pseudospin(r,n,l);
			assert(label[n] == l);
			++l;
		}
	}
	num_loops = l; // set the total number of loops
}

void init_index_shuffle(void)
{
	if (twiceS > 1)
	{
		assert(num_flavour_index_pairs > 0);
		for (size_t i = 0; i < Ns; ++i)
			site_index[i] = i;
	}
	if (twiceS > 2)
	{
		assert(num_flavour_index_pairs > 1);
		auto p = flavour_index_pairs;
		for (size_t b = 1; b < twiceS; ++b)
			for (size_t a = 0; a < b; ++a)
				*p++ = make_pair(a,b);
	}
}

void init_measurements(void)
{
	fill(spin_corr,spin_corr+Ns,0.0);
	fill(spin_corr_collect,spin_corr_collect+Ns,0.0);
}

void symmetrize_fast(size_t bond[D],
               size_t i /*one site*/, 
               size_t a, size_t b /*two local flavours*/)
{
	//cout << "." ;
	assert(a != b);
	const size_t n = twiceS*i+a;
	const size_t m = twiceS*i+b;
	assert(sublat[n] == sublat[m]);
	if (bond[m] == n) return; // already paired
	
	const size_t nn = bond[n];
	assert(bond[nn] == n);
	const size_t mm = bond[m];
	assert(bond[mm] == m);

	const int om = pseudospin[m];
	const int on = pseudospin[n];
	
	// n~~nn       n~~mm  
	//        ==>   
	// m~~mm       m~~nn

	if ((om+on)%2 == 0) // m and n have same pseudospin (randomized loop sublattice structure)
	{
		// flip valence bonds
		bond[n] = mm; 
		bond[mm] = n;
		bond[m] = nn;
		bond[nn] = m;
	}
}

// compute matrix element for strange correlator
inline int spin_half_scorr(size_t n1, size_t n2)
{
	if (label[n1] != label[n2] ) return 0;
	if ((ordering[n1]+ordering[n2]+sublat[n1]+sublat[n2])%2 == 0) return 3;
	return 1;
}

// compute matrix element for (4/3)*(S_i.S_j)
inline int spin_half_corr(size_t n1, size_t n2)
{
	return label[n1] == label[n2] ? 3 : 0;
}
	
double eps_SdotS(size_t i, size_t j)
{
	long int correlation_count = 0; 
	for (size_t a = 0; a < twiceS; ++a)
		for (size_t b = 0; b < twiceS; ++b)
			// normal correlator
			//correlation_count += spin_half_corr(twiceS*i+a,twiceS*j+b);
			// strange correlator
			correlation_count += spin_half_scorr(twiceS*i+a,twiceS*j+b);
	return 0.25*correlation_count;
}

// compute matrix element for (16/3)*(S_i.S_j)(S_k.S_l)
int spin_half_corr(size_t n1, size_t n2, size_t n3, size_t n4)
{
	const int d12 = (label[n1] == label[n2]);
	const int d23 = (label[n2] == label[n3]);
	const int d34 = (label[n3] == label[n4]);		
	int A = 0;
	const int C1 = d12*d23*d34;
	if (C1 != 0)
	{
		const int o1=2*(min(ordering[n1],ordering[n2])/2)+1;
		const int o2=2*(max(ordering[n1],ordering[n2])/2);
		const bool X=(ordering[n3]>=o1 and ordering[n3]<=o2);
		const bool Y=(ordering[n4]>=o1 and ordering[n4]<=o2);
		if (X != Y) A=1;
	}
	const int d13 = (label[n1] == label[n3]);
	const int d24 = (label[n2] == label[n4]);
	const int d14 = (label[n1] == label[n4]);
	const int C2 = d12*d34;
	const int C3 = d13*d24 + d14*d23;
	return -2*(1+2*A)*C1 + 3*C2 + C3;
}

double eps_SdotS(size_t i, size_t j, size_t k, size_t l)
{
	long int correlation_count = 0;
	for (size_t a = 0; a < twiceS; ++a)
		for (size_t b = 0; b < twiceS; ++b)
			for (size_t c = 0; c < twiceS; ++c)
				for (size_t d = 0; d < twiceS; ++d)
			correlation_count += spin_half_corr(twiceS*i+a,twiceS*j+b,twiceS*k+c,twiceS*l+d);
	return (3.0/16.0)*correlation_count;
}

void symmetrize(size_t bond[D])
{
	init_pseudospin(); // use for fast updates
	if (twiceS > 1)
	{
		assert(num_flavour_index_pairs > 0);
		shuffle(site_index,site_index+Ns,rng);
		for (size_t ii = 0; ii < Ns; ++ii)
		{
			const size_t i = site_index[ii];
			
			if (twiceS > 2)
			{
				assert(num_flavour_index_pairs > 1);
				shuffle(flavour_index_pairs,flavour_index_pairs+num_flavour_index_pairs,rng);
				for (size_t f = 0; f < num_flavour_index_pairs; ++f)
				{
					const size_t a = flavour_index_pairs[f].first;
					const size_t b = flavour_index_pairs[f].second;
					symmetrize_fast(bond,i,a,b);
				}
			}
			else
			{
				assert(twiceS == 2);
				assert(num_flavour_index_pairs == 1);
				symmetrize_fast(bond,i,0,1);
			}
		}
	}
	init_labels_fast(); // use for fast updates
}

int main(void)
{
	init_interrupt_handler();
	init_random();
	
	init_config_aklt_square(left);
#ifndef STRANGE
	init_config_aklt_square(right); // Normal correlator for square lattice AKLT state
#else
	init_config_triplet(right); // Strange correlator for square lattice AKLT state
#endif
	
	verify_config(left);
	verify_config(right);
	
	init_index_shuffle();
	init_measurements();
	
	cout << "begin equilibrate" << endl;

	for (size_t e=0; e<1024; ++e) // equilibrate
	{
		symmetrize(left);
		symmetrize(right);
#ifndef NDEBUG
		verify_config(left);
		verify_config(right);
#endif
	}

	cout << "end equilibrate" << endl;

	size_t bin_count = 0;
	size_t num_bins = 0;

	cout << "begin measurement" << endl;

	while(num_bins < 16)
	{
		symmetrize(left);
		symmetrize(right);
#ifndef NDEBUG
		verify_config(left);
		verify_config(right);
#endif
		
		for (size_t i = 0; i < Ns; ++i)
		{
			for (size_t j = 0; j < 1; ++j) // freeze this site
			{
				spin_corr[i] += eps_SdotS(i,j);
			}
		}
		++bin_count;
		if (bin_count%1024 == 0)
		{
			for (size_t i = 0; i < Ns; ++i)
			{
				for (size_t j = 0; j < 1; ++j) // freeze this site
				{
					spin_corr_collect[i] += spin_corr[i]/double(bin_count);
				}
			}
			++num_bins;
			double spin_corr_ave_xx[Lx];
			double spin_corr_ave_xy[Lx];
			double spin_corr_ave_yx[Lx];
			double spin_corr_ave_yy[Lx];
			fill(spin_corr_ave_xx,spin_corr_ave_xx+Lx,0.0);
			fill(spin_corr_ave_xy,spin_corr_ave_xy+Lx,0.0);
			fill(spin_corr_ave_yx,spin_corr_ave_yx+Lx,0.0);
			fill(spin_corr_ave_yy,spin_corr_ave_yy+Lx,0.0);
			for (size_t r = 0; r < Lx; ++r)
			{
				size_t i = r; // along xx direction
				spin_corr_ave_xx[r] += spin_corr_collect[i]/num_bins;
				size_t j = r*(Lx+1)+(Lx-2*r-1)*(r%2); // along xy direction (squares!)
				spin_corr_ave_xy[r] += spin_corr_collect[j]/num_bins;
				size_t k = (r+1)*(Lx-1)-(Lx-2*r-1)*(r%2); // along yx direction (squares!)
				spin_corr_ave_yx[r] += spin_corr_collect[k]/num_bins;
				size_t l = r*Lx+(Lx-1)*(r%2); // along yy direction (squares!)
				spin_corr_ave_yy[r] += spin_corr_collect[l]/num_bins;
			}
			for (size_t r = 0; r < Lx; ++r) // only along one axis
			{
				cout << r
					<< setw(20) << spin_corr_ave_xx[r] / double(1)
					<< setw(20) << spin_corr_ave_xy[r] / double(1)
					<< setw(20) << spin_corr_ave_yx[r] / double(1)
					<< setw(20) << spin_corr_ave_yy[r] / double(1) << endl;
			}
			fill(spin_corr,spin_corr+Ns,0.0);
			bin_count = 0;
			cout << endl;
		}
		
	}
	
	shut_down();
	return 0;	
}
