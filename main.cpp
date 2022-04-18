/*************************************************************************
 SCOPA software:  March, 2016

 Contributors:
 * Andrew P Morris A.P.Morris@liverpool.ac.uk
 * Reedik Magi reedik.magi@ut.ee

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *************************************************************************/

#include <iostream>
#include <algorithm>
#include <string>
#include <vector>
#include <fstream>
#include <cctype> // std::toupper
#include <map>

#include <zlib.h>
#include "CmdLine.h"
#include "global.h"
#include "tools.h"
#include "structures.h"
#include "regression.h"
#include "studenttdistr.h"
#include "chisquaredistr.h"

#include <fstream>
#include <cassert>
#include <stdexcept>
#include <memory>

#include "MissingValue.cpp"
#include "zlib.cpp"
#include "bgen.cpp"

#define LENS 1000000
using namespace TCLAP;
using namespace std;

struct ProbSetter {
	typedef std::vector< std::vector< double > > Data ;
	ProbSetter( Data* result ):
		m_result( result ),
		m_sample_i(0)
	{}

	// Called once allowing us to set storage.
	void initialise( std::size_t number_of_samples, std::size_t number_of_alleles ) {
		m_result->clear() ;
		m_result->resize( number_of_samples ) ;
	}

	// If present with this signature, called once after initialise()
	// to set the minimum and maximum ploidy and numbers of probabilities among samples in the data.
	// This enables us to set up storage for the data ahead of time.
	void set_min_max_ploidy( uint32_t min_ploidy, uint32_t max_ploidy, uint32_t min_entries, uint32_t max_entries ) {
		for( std::size_t i = 0; i < m_result->size(); ++i ) {
			m_result->at( i ).reserve( max_entries ) ;
		}
	}

	// Called once per sample to determine whether we want data for this sample
	bool set_sample( std::size_t i ) {
		m_sample_i = i ;
		// Yes, here we want info for all samples.
		return true ;
	}

	// Called once per sample to set the number of probabilities that are present.
	void set_number_of_entries(
		std::size_t ploidy,
		std::size_t number_of_entries,
		genfile::OrderType order_type,
		genfile::ValueType value_type
	) {
		assert( value_type == genfile::eProbability ) ;
		m_result->at( m_sample_i ).resize( number_of_entries ) ;
		m_entry_i = 0 ;
	}

	// Called once for each genotype (or haplotype) probability per sample.
	void set_value( uint32_t, double value ) {
		m_result->at( m_sample_i ).at( m_entry_i++ ) = value ;
	}

	// Ditto, but called if data is missing for this sample.
	void set_value( uint32_t, genfile::MissingValue value ) {
		// Here we encode missing probabilities with -1
		m_result->at( m_sample_i ).at( m_entry_i++ ) = -1 ;
	}

	// If present with this signature, called once after all data has been set.
	void finalise() {
		// nothing to do in this implementation.
	}

private:
	Data* m_result ;
	std::size_t m_sample_i ;
	std::size_t m_entry_i ;
} ;

// BgenParser is a thin wrapper around the core functions in genfile/bgen/bgen.hpp.
// This class tracks file state and handles passing the right callbacks.
struct BgenParser {

	BgenParser( std::string const& filename ):
		m_filename( filename ),
		m_state( e_NotOpen ),
		m_have_sample_ids( false )
	{
		// Open the stream
		m_stream.reset(
			new std::ifstream( filename, std::ifstream::binary )
		) ;
		if( !*m_stream ) {
			throw std::invalid_argument( filename ) ;
		}
		m_state = e_Open ;

		// Read the offset, header, and sample IDs if present.
		genfile::bgen::read_offset( *m_stream, &m_offset ) ;
		genfile::bgen::read_header_block( *m_stream, &m_context ) ;
		if( m_context.flags & genfile::bgen::e_SampleIdentifiers ) {
			genfile::bgen::read_sample_identifier_block(
				*m_stream, m_context,
				[this]( std::string id ) { m_sample_ids.push_back( id ) ; }
			) ;
			m_have_sample_ids = true ;
		}

		// Jump to the first variant data block.
		m_stream->seekg( m_offset + 4 ) ;

		// We keep track of state (though it's not really needed for this implementation.)
		m_state = e_ReadyForVariant ;
	}

	std::ostream& summarise( std::ostream& o ) const {
		o << "BgenParser: bgen file ("
			<< ( m_context.flags & genfile::bgen::e_Layout ? "v1.2 layout" : "v1.1 layout" )
			<< ", "
			<< ( m_context.flags & genfile::bgen::e_CompressedSNPBlocks ? "compressed" : "uncompressed" ) << ")"
			<< " with "
			<< m_context.number_of_samples << " " << ( m_have_sample_ids ? "named" : "anonymous" ) << " samples and "
			<< m_context.number_of_variants << " variants.\n" ;
		return o ;
	}

	std::size_t number_of_samples() const {
		return m_context.number_of_samples ;
	}

	// Report the sample IDs in the file using the given setter object
	// (If there are no sample IDs in the file, we report a dummy identifier).
	template< typename Setter >
	void get_sample_ids( Setter setter ) {
		if( m_have_sample_ids ) {
			for( std::size_t i = 0; i < m_context.number_of_samples; ++i ) {
				setter( m_sample_ids[i] ) ;
			}
		} else {
			for( std::size_t i = 0; i < m_context.number_of_samples; ++i ) {
				setter( "(unknown_sample_" + std::to_string( i+1 ) + ")" ) ;
			}
		}
	}

	// Attempt to read identifying information about a variant from the bgen file, returning
	// it in the given fields.
	// If this method returns true, data was successfully read, and it should be safe to call read_probs()
	// or ignore_probs().
	// If this method returns false, data was not successfully read indicating the end of the file.
	bool read_variant(
		std::string* chromosome,
		uint32_t* position,
		std::string* rsid,
		std::vector< std::string >* alleles
	) {
		assert( m_state == e_ReadyForVariant ) ;
		std::string SNPID ; // read but ignored in this toy implementation

		if(
			genfile::bgen::read_snp_identifying_data(
				*m_stream, m_context,
				&SNPID, rsid, chromosome, position,
				[&alleles]( std::size_t n ) { alleles->resize( n ) ; },
				[&alleles]( std::size_t i, std::string const& allele ) { alleles->at(i) = allele ; }
			)
		) {
			m_state = e_ReadyForProbs ;
			return true ;
		} else {
			return false ;
		}
	}

	// Read genotype probability data for the SNP just read using read_variant()
	// After calling this method it should be safe to call read_variant() to fetch
	// the next variant from the file.
	void read_probs( std::vector< std::vector< double > >* probs ) {
		assert( m_state == e_ReadyForProbs ) ;
		ProbSetter setter( probs ) ;
		genfile::bgen::read_and_parse_genotype_data_block< ProbSetter >(
			*m_stream,
			m_context,
			setter,
			&m_buffer1,
			&m_buffer2
		) ;
		m_state = e_ReadyForVariant ;
	}

	// Ignore genotype probability data for the SNP just read using read_variant()
	// After calling this method it should be safe to call read_variant()
	// to fetch the next variant from the file.
	void ignore_probs() {
		genfile::bgen::ignore_genotype_data_block( *m_stream, m_context ) ;
		m_state = e_ReadyForVariant ;
	}

private:
	std::string const m_filename ;
	std::unique_ptr< std::istream > m_stream ;

	// bgen::Context object holds information from the header block,
	// including bgen flags
	genfile::bgen::Context m_context ;

	// offset byte from top of bgen file.
	uint32_t m_offset ;

	// We keep track of our state in the file.
	// Not strictly necessary for this implentation but makes it clear that
	// calls must be read_variant() followed by read_probs() (or ignore_probs())
	// repeatedly.
	enum State { e_NotOpen = 0, e_Open = 1, e_ReadyForVariant = 2, e_ReadyForProbs = 3, eComplete = 4 } ;
	State m_state ;

	// If the BGEN file contains samples ids, they will be read here.
	bool m_have_sample_ids ;
	std::vector< std::string > m_sample_ids ;

	// Buffers, these are used as working space by bgen implementation.
	std::vector< genfile::byte_t > m_buffer1, m_buffer2 ;
} ;


double ddabs(double d){if(d<0)return d*-1;return d;}
bool readSampleFile(global & G, ofstream & LOG);
bool readGenoFile(global & G, ofstream & LOG);
bool readExclFile(global & G, ofstream & LOG);
int main (int argc,  char * argv[])
{
    global GLOBAL;

    try
    {
        CmdLine cmd("For more info: http://www.geenivaramu.ee/en/tools/scopa", ' ', GLOBAL.version );
        ValueArg<string> samplefArg("s","sample","This specifies sample file",true,"","string", cmd);
        ValueArg<int> chrArg("", "chr", "This specifies chromosome", false, 0, "int" , cmd);
        ValueArg<string> genofArg("g","gen","This specifies genotype file",true,"","string", cmd);
        ValueArg<string> outfArg("o","out","This specifies output root",true,"","string", cmd);
        ValueArg<string> exclfArg("e","exclusion","This specifies marker exclusion list",false,"","string", cmd);
        ValueArg<string> naArg("","missing_phenotype","This specifies missing data value (default NA)",false,"","string", cmd);
        ValueArg<double> thresholdArg("", "imp_threshold", "Imputation quality threshold (default 0)", false, 0,"double" , cmd);

        MultiArg<string> phenoNamesArg("","pheno_name", "Name of phenotype to use (use this command multiple times i.e. --pheno_name BMI --pheno_name HEIGHT etc.)", true, "string", cmd);

        SwitchArg rmmissingArg("", "remove_missing","Remove sample if any of the phenotype values is missing (default OFF)", cmd);
        SwitchArg printallArg("", "print_all","Print results for all models (default OFF)", cmd);
        SwitchArg printBetas("", "betas","Print effect size and stderr info (default OFF)", cmd);
        SwitchArg printComplexArg("", "print_complex","Print only the model with all phenotypes (default OFF)", cmd);
        SwitchArg printCovarianceArg("", "print_covariance","Print covariance matrix data for the model with all phenotypes (default OFF)", cmd);
        SwitchArg debugArg("", "debug","Debug mode on (default OFF)", cmd);
        cmd.parse(argc,argv);

        GLOBAL.inputSampleFile = samplefArg.getValue();
        GLOBAL.inputGenFile = genofArg.getValue();
        GLOBAL.phenoList = phenoNamesArg.getValue();
        GLOBAL.inputExclFile = exclfArg.getValue();
        if (naArg.getValue() != "")GLOBAL.missingCode = naArg.getValue();
        GLOBAL.removeMissing = rmmissingArg.getValue();
        GLOBAL.printAll = printallArg.getValue();
        GLOBAL.printBetas = printBetas.getValue();
        GLOBAL.printComplex = printComplexArg.getValue();
        GLOBAL.printCovariance = printCovarianceArg.getValue();
        GLOBAL.debugMode = debugArg.getValue();
        GLOBAL.outputRoot = outfArg.getValue();
        GLOBAL.threshold = thresholdArg.getValue();
                GLOBAL.chr = chrArg.getValue();
        if (GLOBAL.phenoList.size()<2)
        {
            cout<< "Less than 2 phenotypes selected for the analysis. Please add additional phenotypes for pleiotropy testing. Exit program!" <<endl;
            exit(1);
        }
        GLOBAL.createOutput();
        ofstream LOG (GLOBAL.outputLog.c_str());
        cout << "###################\n# SCOPA v." << GLOBAL.version << "\n###################\n" << endl;


        LOG << "###################\n# SCOPA v." << GLOBAL.version << "\n###################\n" << endl;
        LOG << "Using following command line options:" << endl;
        LOG << "Genotype file: " << GLOBAL.inputGenFile << endl;
        LOG << "Sample file: " << GLOBAL.inputSampleFile << endl;
        LOG << "Phenotypes: ";
        for (int i=0; i<GLOBAL.phenoList.size(); i++){LOG << GLOBAL.phenoList[i] << " ";}
        LOG << "(" << GLOBAL.phenoList.size() << ")" << endl;
        if (GLOBAL.removeMissing) {LOG << "Remove samples with any missing phenotype data (remove_missing ON)"<<endl;}
        else {LOG << "Use all availabe phenotype data (remove_missing OFF)"<<endl;}
        LOG << "Output result file: " << GLOBAL.outputResult << endl;
        LOG << "Output log file: " << GLOBAL.outputLog << endl;
        LOG << "Output betas file: " << GLOBAL.outputBetas << endl;

        if (GLOBAL.printCovariance && GLOBAL.printComplex)
        {
            LOG << "Print covariance matrix (print_covariance ON)"<<endl;
        }
        else if (GLOBAL.printCovariance && !GLOBAL.printComplex)
        {
            cout << "printCovariance option is only available, if printComplex is used. Exit program.";
            exit(1);
        }
        if (GLOBAL.printAll) {LOG << "Print all possible models (print_all ON)"<<endl;}
        else if (GLOBAL.printComplex){LOG << "Print only the model with all phenotypes (print_complex ON)"<<endl;}
        else {LOG << "Print only best model (print_all OFF)"<<endl;}
        if (GLOBAL.debugMode) {LOG << "DEBUG MODE ON"<<endl;}
        if (GLOBAL.printBetas){LOG << "Creating file for saving all beta and stderr values for all phenotypes in all selected models" << endl;}

        if (GLOBAL.threshold<0 || GLOBAL.threshold>1)
        {
            cout << "Imputation quality threshold out of range. Must be between 0  and 1. Exit program.";
            exit(1);
        }

        LOG << "Imputation quality threshold: " << GLOBAL.threshold << endl;

        cout << "Reading sample file..." << endl;
        readSampleFile(GLOBAL, LOG);
        if (GLOBAL.inputExclFile != "")
        {
            cout << "Reading exclusion list file..." << endl;
            readExclFile(GLOBAL, LOG);
        }
        cout << "Reading genotype file..." << endl;
        readGenoFile(GLOBAL, LOG);

        LOG << "Analysis finished" <<endl;




    } catch (ArgException &e)  // catch any exceptions
    { cerr << "error: " << e.error() << " for arg " << e.argId() << endl; }


    return 0;
}



bool
readExclFile(global & G, ofstream & LOG)
{
    int lineNr = 0;
    ifstream F (G.inputExclFile.c_str());
    if (F.is_open())
    {
       	while (! F.eof() )
        {
	        string line;
        	vector<string> tokens;
        	getline (F,line);
			string currentmarker = "";
			int n = Tokenize(string(line), tokens, " ");		//tabulating file by space
			if (n>0)
			{
                G.exclusionList[string(tokens[0])]= true;
			}
            lineNr++;
		}
    }
    else
    {cout << "Cannot exclusion list file. Exit program!" << endl;exit(1);}
    if (G.debugMode)cout << "Exclusion list contained: " << lineNr << " markers" << endl;
    LOG << "Exclusion list contained: " << lineNr << " markers" << endl;
    return true;
}


bool
readSampleFile(global & G, ofstream & LOG)
{
    int lineNr = 0;
    ifstream F (G.inputSampleFile.c_str());
    if (F.is_open())
    {
        vector <string> columnNames;
        vector <int> phenoColumns;
        vector <string> sampleNames;
        while (! F.eof() )
        {
            string line;
            vector<string> tokens;
            getline (F,line);
            string currentmarker = "";
            int n = Tokenize(string(line), tokens, " ");            //tabulating file by space
            if (lineNr==0)
            {
                if (n>2)
                {
                    for (int i = 0; i < n; i++)
                    {
                        for (int j = 0; j<G.phenoList.size(); j++)
                        {
                            if (uc(tokens[i]) == uc(G.phenoList[j]))
                            {
                                phenoColumns.push_back(i);
                                columnNames.push_back(tokens[i]);
                            }
                        }
                    }

                }
                if (G.phenoList.size()!=phenoColumns.size())
                {
                    cout << "One or more phenotypes cannot be found from the sample file. Exit program!\nExpected: ";
                    for (int j = 0; j<G.phenoList.size(); j++){cout << G.phenoList[j] << " ";}
                    cout << "\nFound: ";
                    for (int j = 0; j<phenoColumns.size(); j++){cout << columnNames[j] << " ";}
                    cout << endl;
                    exit(1);
                }
                else
                {
                    G.phenoList = columnNames;
                }
            }
            if (lineNr>=2)   // read data
            {
                if (n>2)
                {
                    vector <double> _phenos;
                    string _name = tokens[0];
                    bool isOK=true;
                    bool _hasMissing = false;
                    for (int i = 0; i < phenoColumns.size(); i++)
                    {
                        if (tokens[phenoColumns[i]] == G.missingCode){_phenos.push_back(-9999);_hasMissing=true;}
                        else _phenos.push_back(atof(tokens[phenoColumns[i]].c_str()));
                    }
                    if (G.removeMissing==true && _hasMissing==true)
                    {
                        isOK=false;
                        for (int i = 0; i < phenoColumns.size(); i++){_phenos[i]=-9999;}
                    }
                    sample _S;
                    _S._name=_name;
                    _S._phenos = _phenos;
                    _S.isOK = isOK;
                    G.samples.push_back(_S);
                    /*if (G.debugMode) cout << _name << ":\t";
                    for (int j = 0; j<G.phenoList.size(); j++)
                    {
                      if (G.debugMode)cout << _phenos[j] << "\t";
                    }
                    if (G.debugMode)cout << endl;
                    */
                }
            }
            lineNr++;
        }
				// cout << "sample size: " << G.samples.size() << endl;
    }
    else
    {cout << "Cannot read sample file. Exit program!" << endl;exit(1);}
    if (G.debugMode)cout << "Altogether: " << G.phenoList.size() << " phenos" << endl;
    if (G.debugMode)cout << "Altogether: " << G.samples.size() << " samples" << endl;
    LOG << "Sample file contained all: " << G.phenoList.size() << " phenotypes" << endl;
    LOG << "Sample file contained: " << G.samples.size() << " samples" << endl;
    return true;
}

bool
readGenoFile(global & G, ofstream & LOG)
{
    ofstream OUT (G.outputResult.c_str());
    ofstream BETAS (G.outputBetas.c_str());
    if (!G.printCovariance)
    OUT << "Chromosome\tPosition\tMarkerName\tEffectAllele\tOtherAllele\tInfoScore\tHWE\tMAF\tN\tAA\tAB\tBB\tPhenotypeCount\tMask\tLogLikelihood\tnullLogLikelihood\tLikelihoodRatio\tP-value\tBIC\tBICnull\tModel\tsortedModel\n";
    else
    {
        OUT << "Chromosome\tPosition\tMarkerName\tEffectAllele\tOtherAllele\tInfoScore\tHWE\tMAF\tN\tAA\tAB\tBB\tPhenotypeCount\tMask\tLogLikelihood\tnullLogLikelihood\tLikelihoodRatio\tP-value\tBIC\tBICnull\tModel\tsortedModel";
        for (int i = 0; i < G.phenoList.size();i++) OUT << "\tbeta_" << i+1 << "\tse_" << i+1;
        for (int i = 0; i < G.phenoList.size();i++)
        {
            for (int j = i; j < G.phenoList.size();j++)
                OUT << "\tcov_" << i+1 << "_" << j+1;
        }
        OUT << endl;
    }
    if (G.printBetas)BETAS << "MarkerName\tEffectAllele\tOtherAllele\tN\tModel\tModel_member\tbeta\tse\n";

		// Try reading BGEN file
		if (G.inputGenFile.substr(G.inputGenFile.length()-4)=="bgen")
		{
			try{
				string const filename = G.inputGenFile;
				BgenParser bgenParser(filename) ;

				// To be removed, for debugging
				// Summarise information about the bgen file
				bgenParser.summarise(cerr) ;

				// Print out all sample ids, if given
				// bgenParser.get_sample_ids(
				// 	[](string const& id ) {cout << "\t" << id ; }
				// ) ;

				// To store what's given by the function read_variant()
				string chromosome ;
				uint32_t position ;
				string rsid ;
				vector<string> alleles ;
				vector<vector<double>> probs ;

				// VARIANT
		    while( bgenParser.read_variant( &chromosome, &position, &rsid, &alleles )){

					// CHROMOSOME
					int chr; // To be printed out in debug mode
					if (chromosome=="MT") chr=26;
					else if (chromosome=="XY") chr=25;
					else if (chromosome=="Y" or chromosome == "0Y") chr=24; // Have to specify "0X" and "0Y" for some syntax variation
					else if (chromosome=="X" or chromosome == "0X") chr=23;
					else chr = atoi(chromosome.c_str());

					if (G.chr)chr=G.chr;
					if (G.debugMode) cout << "Chromosome id: " << chr;

					// POSITION, MARKER
					std::cout << chr << ' ' << rsid << ' ' << position << ' '; // To be removed
					if (G.debugMode) cout << "Pos: " << position << "\nmarker:" << rsid;

					// EFFECT & NON-EFFECT ALLELES
					// To be removed
					assert( alleles.size() > 0 ) ;
					std::cout << alleles[0] << ' ' ;
					for( std::size_t i = 1; i < alleles.size(); ++i ) {
						std::cout << alleles[i] ;
					}

					// Debug mode
					string effectAllele;
					string nonEffectAllele;
					if (alleles.size() == 2){
						effectAllele = alleles[1];
						nonEffectAllele = alleles[0];
					}
					if (G.debugMode) cout << "\nea/nea:" << effectAllele <<"/" << nonEffectAllele<<"\n";

					// bgenParser.ignore_probs();

					// PROBABILITIES
					// Initialise variables
					vector <double> gen;
					vector <double> gen2;
					int j = 0;
					double aa=0; double aA=0; double AA=0;
					double callrate=0; double ok_gen=0; double not_ok_gen=0;
					double bestModel = 1e200;
					string bestModelString = "";
					string bestBetasString = "";
					double fijeij = 0; //for infoscore
					double infoscore = 1;
					double gen2_pb_val = 0; double gen_pb_val = 0;
					double eij = 0; double fij = 0;

					// CALL RATE
					// Read probability
					bgenParser.read_probs( &probs );
					for( size_t i = 0; i < probs.size(); ++i ) {

						// For each sample/individual, reset some variables
						cout << ' ' ; // To be removed
						gen2_pb_val = 0;
						gen_pb_val = 0;
						eij = 0;
						fij = 0;

						for( size_t j = 0; j < probs[i].size(); ++j ) {
							if( probs[i][j] == -1 ) {
								cout << "." ; // Probability data unavailable
							}
							else{
								cout << probs[i][j] << " "; // Probability data available

								switch (j) {
									case 0:
									aa+=probs[i][j];
									gen_pb_val+=2*(probs[i][j]);
									case 1:
									aA+=probs[i][j];
									gen_pb_val+=(probs[i][j]);
									gen2_pb_val+=(probs[i][j]);
									eij+=probs[i][j];
									fij+=probs[i][j];
									case 2:
									AA+=probs[i][j];
									gen2_pb_val+=2*(probs[i][j]);
									eij+=2*(probs[i][j]);
									fij+=4*(probs[i][j]);
								}
							}
						}

						ok_gen++;
						gen2.push_back(gen2_pb_val);
						gen.push_back(gen_pb_val);
						fijeij+= fij - (eij*eij);
						j++; // what is this for?
					}

					if (ok_gen+not_ok_gen!=G.samples.size())
					{
							cout << "The number of samples in genotype file (" << ok_gen+not_ok_gen << ") does not match the number of samples in sample file (" << G.samples.size() << "). Exit program!" << endl;
							exit (1);
					}
					callrate = ok_gen/(ok_gen+not_ok_gen);
					if (G.debugMode)cout<<"Callrate: "<< callrate <<endl;

					// Print to check, to be removed
					cout << "aa: " << aa << " aA: " << aA << " AA: " << AA;
					cout << "\n";
					cout << "eij: " << eij << " fij: " << fij;
					cout << "\n";
					cout << "fijeij: " << fijeij;
					cout << '\n';
				}

				return 0;
			}
			catch( genfile::bgen::BGenError const& e ) {
				std::cerr << "!! Uh-oh, error parsing bgen file.\n" ;
				return -1 ;
			}
		}


    if (G.inputGenFile.substr(G.inputGenFile.length()-2)=="gz")
    {
        //      	ifstream F (G.inputGenFile.c_str());
        gzFile F =gzopen(G.inputGenFile.c_str(),"r");
        char *buffer = new char[LENS];
        while(0!=gzgets(F,buffer,LENS))
        {
                string line;
                vector<string> tokens;

                string currentmarker = "";
                int n = Tokenize(buffer, tokens, " ");            //tabulating file by space
                if (n>2 && !G.exclusionList[string(tokens[1])])
                {
                    int chr;
                    if (uc(tokens[0])=="MT") chr=26;
                    else if (uc(tokens[0])=="XY") chr=25;
                    else if (uc(tokens[0])=="Y") chr=24;
                    else if (uc(tokens[0])=="X") chr=23;
                    else chr = atoi(tokens[0].c_str());
                    if (G.chr)chr=G.chr;
                    if (G.debugMode) cout << "Chromosome id: " << chr;
                    int pos = atoi(tokens[2].c_str());
                    string markerName = string(tokens[1]);
                    string effectAllele = tokens[4];
                    string nonEffectAllele = tokens[3];
                    if (G.debugMode) cout << "Pos: " << pos << "\nmarker:" << markerName << "\nea/nea:" << effectAllele <<"/" << nonEffectAllele<<"\n";

                    vector <double> gen;
                    vector <double> gen2;
                    int j = 0;
                    double aa=0; double aA=0; double AA=0;
                    double callrate=0; double ok_gen=0; double not_ok_gen=0;
                    double bestModel = 1e200;
                    string bestModelString = "";
                    string bestBetasString = "";
                    double fijeij = 0; //for infoscore
                    double infoscore = 1;
                    for (int i = 5; i < n-1; i+=3)
                    {
                        //granvil copypaste
                        aa+=atof(tokens[i].c_str());
                        aA+=atof(tokens[i+1].c_str());
                        AA+=atof(tokens[i+2].c_str());
                        gen2.push_back(((atof(tokens[i+2].c_str())*2) + atof(tokens[i+1].c_str())));
                        gen.push_back(((atof(tokens[i].c_str())*2) + atof(tokens[i+1].c_str())));
                        ok_gen++;
                        double eij=(2*atof(tokens[i+2].c_str())) + atof(tokens[i+1].c_str());
                        double fij = (4*atof(tokens[i+2].c_str())) + atof(tokens[i+1].c_str());
                        fijeij+= fij - (eij*eij);



                        j++;
                    }
                                        if (ok_gen+not_ok_gen!=G.samples.size())
                    {
                        cout << "The number of samples in genotype file (" << ok_gen+not_ok_gen << ") does not match the number of samples in sample file (" << G.samples.size() << "). Exit program!" << endl;
                        exit (1);
                    }
                    callrate = ok_gen/(ok_gen+not_ok_gen);
                    if (G.debugMode)cout<<"Callrate: "<< callrate <<endl;
                    double maf = 0;
                    if (aa+aA+AA>0)
                    {
                        maf = ((2*aa)+aA)/(2*(aa+aA+AA));

                    }
                    if (maf==0){infoscore=1;}
                    //info score calculation according to the measure in snptest:
                    //https://mathgen.stats.ox.ac.uk/genetics_software/snptest/snptest.v2.pdf
                    else
                    {
                        infoscore = 1-(fijeij/(2*(aa+aA+AA)*maf*(1-maf)));
                    }

                    bool firstIsMajorAllele = false;
                    if (maf>0.5)
                    {
                        maf=1-maf;
                        firstIsMajorAllele=true;
                        effectAllele = tokens[3];
                        nonEffectAllele = tokens[4];
                    }
                    if (G.debugMode)cout<<"MAF: "<< maf <<endl;
                    if (maf>0) //marker is ok - Sham'st thou to show thy dangerous brow by night, When evils are most free?
                    {
                        infoscore = 1-(fijeij/(2*(aa+aA+AA)*maf*(1-maf)));
                        if (infoscore >= G.threshold)
                        {
                            unsigned int _rows = (n-5)/3;
                            unsigned int _cols = (unsigned int) G.phenoList.size()+1;
                            if (G.debugMode)cout << "Mainmatrix x: " << _rows << " y: " <<  _cols << endl;
                            matrixD _mainmatrix(_rows, _cols); //[0]-Y(0,1,2,3); [1]-.. X1,X2,X3...
                            int curind=0;
                            for (int i = 5; i < n-1; i+=3)
                            {
                                if (firstIsMajorAllele)
                                {
                                    _mainmatrix.put(curind,0,((atof(tokens[i].c_str())*2) + atof(tokens[i+1].c_str())));
                                }
                                else
                                {
                                    _mainmatrix.put(curind,0,((atof(tokens[i+2].c_str())*2) + atof(tokens[i+1].c_str())));

                                }
                                for (int j = 0; j<G.phenoList.size(); j++)
                                {
                                    string xxy = G.samples[curind]._name;
                               //     int xxyy = j;
                                 //   double xxx = G.samples[curind]._phenos[j];

                                    _mainmatrix.put(curind,j+1,G.samples[curind]._phenos[j]);
                                }
                                curind++;
                            }
                            //We have main table now - lets run through all possible combinations of phenotypes
                            //if (G.debugMode)_mainmatrix.print();
                            int _testcount = pow(2,G.phenoList.size());
                            for (int test=_testcount-1; test>=1; test--)
                            {
                                int _phenoCount = 0;
                                int _sampleCount = 0;
                                vector<bool> phenoMask = phenoMasker(test,(int)G.phenoList.size());

                                if (G.debugMode)cout << "Current mask: ";
                                for (int i = 0; i < phenoMask.size(); i++)
                                {
                                    if (phenoMask[i])
                                    {
                                        _phenoCount++;
                                        if (G.debugMode)cout<<"1";
                                    }
                                    else{if (G.debugMode)cout<<"0";}
                                }
                                if (G.debugMode)cout << endl;


                                for (int i = 0; i<_mainmatrix.getRows(); i++)
                                {
                                    bool indISOK = true;
                                    if (G.debugMode){cout << G.samples[i]._name;}
                                    if (_mainmatrix.get(i,0)==-9999){indISOK=false;}

                                    for (int j = 1; j<_mainmatrix.getCols(); j++)
                                    {

                                        //phenotypes
                                        if (phenoMask[j-1]){
                                            //if (i==0)_phenoCount++;
                                            if (_mainmatrix.get(i,j)==-9999){indISOK=false;}
                                        }
                                        if (G.debugMode){cout << " " << _mainmatrix.get(i,j);}

                                    }
                                    if (G.debugMode){cout << endl;}

                                    if (indISOK)_sampleCount++;
                                }
                                if (G.debugMode)cout << "test: " << test << " indcount: " << _sampleCount << " phenocount: " <<   _phenoCount << endl;
                                arrayD *Y = new arrayD ( _sampleCount );
                                matrixD *X = new matrixD ( _sampleCount, _phenoCount+1);
                                arrayD *W = new arrayD ( _sampleCount );
                                int curInd = 0;
                                vector <double> COPYpheno;
                                for (int i = 0; i<_mainmatrix.getRows(); i++)
                                {
                                    bool indISOK = true;
                                    if (_mainmatrix.get(i,0)==-9999){indISOK=false;}
                                    for (int j = 1; j<_mainmatrix.getCols(); j++)
                                    {

                                        if (phenoMask[j-1]){
                                            //if (i==0)_phenoCount++;
                                            if (_mainmatrix.get(i,j)==-9999){indISOK=false;}
                                        }
                                    }
                                    if (indISOK)
                                    {
                                        int z=1;

                                        //genotype
                                        Y->put(curInd, _mainmatrix.get(i,0));
                                        COPYpheno.push_back(_mainmatrix.get(i,0));

                                        W->put(curInd,1);
                                        X->put(curInd, 0,  1.0);

                                        for (int j = 1; j<_mainmatrix.getCols(); j++)
                                        {
                                            if (phenoMask[j-1])
                                            {
                                                X->put(curInd, z, _mainmatrix.get(i,j));
                                                z++;
                                            }
                                        }
                                        curInd++;
                                    }
                                }
                                //if (G.debugMode)cout << "Print Y" << endl;
                                //if (G.debugMode)Y->print();
                                //if (G.debugMode)cout << "Print X" << endl;
                                //if (G.debugMode)X->print();
                                //if (G.debugMode)
                                //{
                                //    for (int xxx=0; xxx<COPYpheno.size(); xxx++)cout<< COPYpheno[xxx];
                                //    cout <<endl;
                                //}

                                lr LR;
                                try {
                                    if (LR.lr_w(Y, X, W) && test!=_testcount)
                                    {
                                        matrixD * C = LR.covariance;
                                        C->InvertS();
                                        double x;
                                        for (int i=0; i<(unsigned int) _phenoCount;i++)
                                            for (int j=0; j<(unsigned int) _phenoCount;j++)
                                                x+=LR.Cstat->get(i)*LR.Cstat->get(j)*C->get(i,j);
                                        double testLogLikelihood = -LR.getLnLk(Y, X, W)/2;
                                        double nullLogLikelihood = -LR.nullLikelihood(COPYpheno)/2;
                                        if (G.debugMode){cout << "Linear regression done with model: " << test << endl;}
                                        for (int k = 0; k<=_phenoCount;k++)
                                        {
                                            double _pPheno = (1-studenttdistribution((_sampleCount*2) - _phenoCount,ddabs(LR.Cstat->get(k)/LR.SECstat->get(k))))*2;
                                            if (G.debugMode){cout << k << "\t" << (_sampleCount*2) << "\t" << _phenoCount << "\t" << LR.Cstat->get(k) <<  "\t" << LR.SECstat->get(k) << "\t" << _pPheno << endl;}

                                        }
                                        double likelihoodRatio = 2 * (testLogLikelihood - nullLogLikelihood);
                                        double _BIC = (-2 * testLogLikelihood) + ((_phenoCount+1) * log(_sampleCount));
                                        double _BICnull = (-2 *nullLogLikelihood) + (log(_sampleCount));
                                        //double _pModel = gsl_cdf_chisq_Q( likelihoodRatio , _phenoCount); // Overall P value - need to figure out the number of df!!!!!
                                        double _pModel;
                                        if (ddabs(likelihoodRatio)>0){_pModel = 1-chisquaredistribution(_phenoCount,ddabs(likelihoodRatio));}
                                        else {_pModel = NAN;}

                                        if (G.debugMode){cout << "Likelihood: " << testLogLikelihood << endl;}
                                        if (G.debugMode){cout << "nullLikelihood: " << nullLogLikelihood << endl;}
                                        if (G.debugMode){cout << "Model likelihood ratio: " << likelihoodRatio << endl;}
                                        if (G.debugMode){cout << "Model p: " << _pModel << endl;}

                                        /*
                                         cout << chr << "\t"<<  pos << "\t" <<  markerName << "\t" << effectAllele <<"\t" << nonEffectAllele<< "\t" << _phenoCount << "\t";
                                         for (int i = 0; i < phenoMask.size(); i++) {if (phenoMask[i]){cout<<"1";}else{cout<<"0";}}
                                         cout << "\t" << LR.getLnLk(Y, X, W) <<  "\t" << -2 * nullLikelihood(COPYpheno) <<
                                         "\t" << likelihoodRatio << "\t" << _pModel << "\t" << _BIC << "\t" << _BICnull << "\t";
                                         bool phenostart = true;
                                         for (int i = 0; i < phenoMask.size(); i++) {if (phenoMask[i] && phenostart){cout<<G.phenoList[i];phenostart=false;}else if(phenoMask[i]){cout << "+" <<G.phenoList[i];}}
                                         cout << endl;
                                         */
                                        std::stringstream line, line2;
                                        if (_BIC < bestModel)
                                        {

                                            line << chr << "\t"<<  pos << "\t" <<  markerName << "\t" << effectAllele <<"\t" << nonEffectAllele<< "\t" << infoscore << "\t" << HWE(aa,aA,AA) << "\t" <<maf << "\t" << _sampleCount << "\t";
                                            if (firstIsMajorAllele){line << AA << "\t" << aA << "\t" << aa;}
                                            else {line << aa << "\t" << aA << "\t" << AA;}
                                            line << "\t" << _phenoCount << "\t";
                                            for (int i = 0; i < phenoMask.size(); i++) {if (phenoMask[i]){line<<"1";}else{line<<"0";}}
                                            line << "\t" << testLogLikelihood<<  "\t" <<  nullLogLikelihood <<
                                            "\t" << likelihoodRatio << "\t" << _pModel << "\t" << _BIC << "\t" << _BICnull << "\t";
                                            //bool phenostart = true;
                                            //for (int i = 0; i < phenoMask.size(); i++) {if (phenoMask[i] && phenostart){line<<G.phenoList[i];phenostart=false;}else if(phenoMask[i]){line << "+" <<G.phenoList[i];}}
                                            vector<string> phenosorter;
                                            for (int i = 0; i < phenoMask.size(); i++){if(phenoMask[i])phenosorter.push_back(G.phenoList[i]);}
                                            bool phenostart = true;
                                            for (int i = 0; i < phenosorter.size(); i++) {if (phenostart){line<<phenosorter[i];phenostart=false;}else{line << "+" <<phenosorter[i];}}
                                            line << "\t";
                                            sort(phenosorter.begin(),phenosorter.end());
                                            phenostart = true;
                                            for (int i = 0; i < phenosorter.size(); i++) {if (phenostart){line<<phenosorter[i];phenostart=false;}else{line << "+" <<phenosorter[i];}}

                                            if (G.printCovariance)
                                            {
                                                for (int i = 1; i < LR.Cstat->size() ;i++) line << "\t" << LR.Cstat->get(i) << "\t" <<LR.SECstat->get(i);
                                                for (int i = 1; i < LR.Cstat->size();i++)
                                                {
                                                    for (int j = i; j < LR.Cstat->size();j++)
                                                        line << "\t" << LR.covariance->get(i,j);
                                                }
                                            }
                                            line << endl;
                                            bestModelString = line.str();

                                            int k=1;
                                            for (int i = 0; i < phenoMask.size(); i++)
                                            {
                                                if (phenoMask[i])
                                                {
                                                    line2 << markerName << "\t" << effectAllele <<"\t" << nonEffectAllele<< "\t" <<  _sampleCount << "\t";
                                                    bool phenostart = true;
                                                    for (int j = 0; j < phenoMask.size(); j++)
                                                    {
                                                        if (phenoMask[j] && phenostart)
                                                        {
                                                            line2<<G.phenoList[j];
                                                            phenostart=false;
                                                        }
                                                        else if(phenoMask[j])
                                                        {
                                                            line2 << "+" <<G.phenoList[j];
                                                        }
                                                    }
                                                    line2  << "\t" << G.phenoList[i] << "\t" << LR.Cstat->get(k) << "\t" << LR.SECstat->get(k) << endl;
                                                    k++;
                                                }
                                            }
                                            bestBetasString = line2.str();
                                            bestModel = _BIC;
                                        }
                                        if (G.printAll || G.printComplex)
                                        {

                                            OUT << chr << "\t"<<  pos << "\t" <<  markerName << "\t" << effectAllele <<"\t" << nonEffectAllele<< "\t" << infoscore << "\t" << HWE(aa,aA,AA) << "\t" <<maf << "\t" << _sampleCount << "\t";
                                            if (firstIsMajorAllele){OUT << AA << "\t" << aA << "\t" << aa;}
                                            else {OUT << aa << "\t" << aA << "\t" << AA;}

                                            OUT << "\t" << _phenoCount << "\t";
                                            for (int i = 0; i < phenoMask.size(); i++) {if (phenoMask[i]){OUT<<"1";}else{OUT<<"0";}}
                                            OUT << "\t" << testLogLikelihood <<  "\t" <<  nullLogLikelihood <<
                                            "\t" << likelihoodRatio << "\t" << _pModel << "\t" << _BIC << "\t" << _BICnull << "\t";
                                            //bool phenostart = true;
                                            //for (int i = 0; i < phenoMask.size(); i++) {if (phenoMask[i] && phenostart){OUT<<G.phenoList[i];phenostart=false;}else if(phenoMask[i]){OUT << "+" <<G.phenoList[i];}}
                                            vector<string> phenosorter;
                                            for (int i = 0; i < phenoMask.size(); i++){if(phenoMask[i])phenosorter.push_back(G.phenoList[i]);}
                                            bool phenostart = true;
                                            for (int i = 0; i < phenosorter.size(); i++) {if (phenostart){OUT<<phenosorter[i];phenostart=false;}else{OUT << "+" <<phenosorter[i];}}
                                            OUT << "\t";
                                            sort(phenosorter.begin(),phenosorter.end());
                                            phenostart = true;

                                            for (int i = 0; i < phenosorter.size(); i++) {if (phenostart){OUT<<phenosorter[i];phenostart=false;}else{OUT << "+" <<phenosorter[i];}}

                                            if (G.printCovariance)
                                            {
                                                for (int i = 1; i < LR.Cstat->size() ;i++) OUT << "\t" << LR.Cstat->get(i) << "\t" <<LR.SECstat->get(i);
                                                for (int i = 1; i < LR.Cstat->size();i++)
                                                {
                                                    for (int j = i; j < LR.Cstat->size();j++)
                                                        OUT << "\t" << LR.covariance->get(i,j);
                                                }
                                            }
                                            OUT << endl;
                                            int k=1;
                                            for (int i = 0; i < phenoMask.size(); i++)
                                            {
                                                if (phenoMask[i])
                                                {
                                                    if (G.printBetas) BETAS << markerName << "\t" << effectAllele <<"\t" << nonEffectAllele<< "\t" <<  _sampleCount << "\t";
                                                    bool phenostart = true;
                                                    for (int j = 0; j < phenoMask.size(); j++)
                                                    {
                                                        if (phenoMask[j] && phenostart)
                                                        {
                                                            if (G.printBetas) BETAS<<G.phenoList[j];
                                                            phenostart=false;
                                                        }
                                                        else if(phenoMask[j])
                                                        {
                                                            if (G.printBetas) BETAS << "+" <<G.phenoList[j];
                                                        }
                                                    }
                                                    if (G.printBetas) BETAS  << "\t" << G.phenoList[i] << "\t" << LR.Cstat->get(k) << "\t" << LR.SECstat->get(k) << endl;
                                                    k++;
                                                }
                                            }
                                        }

                                        /*else if ((LR.getLnLk(Y, X, W) <= (-2 * nullLikelihood(COPYpheno))) && ((-2 * nullLikelihood(COPYpheno)) > bestModel))
                                         {
                                         line << chr << "\t"<<  pos << "\t" <<  markerName << "\t" << effectAllele <<"\t" << nonEffectAllele<< "\t" << 0 << "\t";
                                         for (int i = 0; i < phenoMask.size(); i++) {if (phenoMask[i]){line<<"0";}else{line<<"0";}}
                                         line << "\t" <<-2 * nullLikelihood(COPYpheno)  <<  "\t" << -2 * nullLikelihood(COPYpheno) <<
                                         "\t" << 1 << "\t" << 1 << "\t" << _BICnull << "\t" << _BICnull << "\tnull\n";
                                         bestModelString = line.str();
                                         bestModel = (-2 * nullLikelihood(COPYpheno));
                                         }*/
                                        // cout << setprecision(4) << test << " " << Pmodel << " "
                                        // << LR.Cstat->get(3) << " " << LR.SECstat->get(3) << " " << P1 << " "
                                        // << LR.Cstat->get(1) << " " << LR.SECstat->get(2) << " " << P2 << " "
                                        // << LR.Cstat->get(2) << " " << LR.SECstat->get(2) << " " << P3 << "\n" ;
                                    }
                                    else
                                    {
                                        if (test!=_testcount) LOG << "Collinearity problem with model: " << chr << "\t"<<  pos << "\t" <<  markerName << "\t" << test << " " << _testcount << " ";
                                        bool phenostart = true;
                                        for (int i = 0; i < phenoMask.size(); i++) {if (phenoMask[i] && phenostart){LOG<<G.phenoList[i];phenostart=false;}else if(phenoMask[i]){LOG << "+" <<G.phenoList[i];}}
                                        LOG << endl;
                                    }

                                } catch (ArgException &e)  // catch any exceptions
                                { cerr << "error: " << e.error() << " for linear model " << e.argId() << endl; }
                                delete X;
                                delete Y;
                                delete W;
                                if (G.printComplex){break;}
                            }

                            if (!G.printAll && !G.printComplex)
                            {
                                if (G.printBetas)BETAS << bestBetasString;
                                OUT << bestModelString;
                            }
                        } //infoscore end i think

                    } //maf > 0 end (i think)
            }
        }

    }
    else
    {
        ifstream F (G.inputGenFile.c_str());
        if (F.is_open())
        {
            while (! F.eof() )
            {
                string line;
                vector<string> tokens;
                getline (F,line);
                string currentmarker = "";
                int n = Tokenize(string(line), tokens, " ");            //tabulating file by space
                if (n>2 && !G.exclusionList[string(tokens[1])])
                {
                    int chr;
                    if (uc(tokens[0])=="MT") chr=26;
                    else if (uc(tokens[0])=="XY") chr=25;
                    else if (uc(tokens[0])=="Y") chr=24;
                    else if (uc(tokens[0])=="X") chr=23;
                    else chr = atoi(tokens[0].c_str());
                    if (G.chr)chr=G.chr;
                    if (G.debugMode) cout << "Chromosome id: " << chr;
                    int pos = atoi(tokens[2].c_str());
                    string markerName = string(tokens[1]);
                    char effectAllele = tokens[4][0];
                    char nonEffectAllele = tokens[3][0];
                    if (G.debugMode) cout << "Pos: " << pos << "\nmarker:" << markerName << "\nea/nea:" << effectAllele <<"/" << nonEffectAllele<<"\n";

                    vector <double> gen, gen2;
                    int j = 0;
                    double aa=0; double aA=0; double AA=0;
                    double callrate=0; double ok_gen=0; double not_ok_gen=0;
                    double bestModel = 1e200;
                    string bestModelString = "";
                    string bestBetasString = "";
                    double fijeij = 0; //for infoscore
                    double infoscore = 1;
                    for (int i = 5; i < n-1; i+=3)
                    {
                        //granvil copypaste
                        aa+=atof(tokens[i].c_str());
                        aA+=atof(tokens[i+1].c_str());
                        AA+=atof(tokens[i+2].c_str());
                        gen2.push_back(((atof(tokens[i+2].c_str())*2) + atof(tokens[i+1].c_str())));
                        gen.push_back(((atof(tokens[i].c_str())*2) + atof(tokens[i+1].c_str())));
                        ok_gen++;
                        double eij=(2*atof(tokens[i+2].c_str())) + atof(tokens[i+1].c_str());
                        double fij = (4*atof(tokens[i+2].c_str())) + atof(tokens[i+1].c_str());
                        fijeij+= fij - (eij*eij);



                        j++;
                    }
                    if (ok_gen+not_ok_gen!=G.samples.size())
                    {
                        cout << "The number of samples in genotype file (" << ok_gen+not_ok_gen << ") does not match the number of samples in sample file (" << G.samples.size() << "). Exit program!" << endl;
                        exit (1);
                    }
                    callrate = ok_gen/(ok_gen+not_ok_gen);
                    if (G.debugMode)cout<<"Callrate: "<< callrate <<endl;
                    double maf = 0;
                    if (aa+aA+AA>0)
                    {
                        maf = ((2*aa)+aA)/(2*(aa+aA+AA));

                    }
                    bool firstIsMajorAllele = false;
                    if (maf>0.5)
                    {
                        maf=1-maf;
                        firstIsMajorAllele=true;
                        effectAllele = tokens[3][0];
                        nonEffectAllele = tokens[4][0];
                    }
                    if (G.debugMode)cout<<"MAF: "<< maf <<endl;
                    if (maf>0) //marker is ok - Sham'st thou to show thy dangerous brow by night, When evils are most free?
                    {
                        infoscore = 1-(fijeij/(2*(aa+aA+AA)*maf*(1-maf)));
                        if (infoscore>=G.threshold)
                        {
                            unsigned int _rows = (n-5)/3;
                            unsigned int _cols = (unsigned int) G.phenoList.size()+1;
                            if (G.debugMode)cout << "Mainmatrix x: " << _rows << " y: " <<  _cols << endl;
                            matrixD _mainmatrix(_rows, _cols); //[0]-Y(0,1,2,3); [1]-.. X1,X2,X3...
                            int curind=0;
                            for (int i = 5; i < n-1; i+=3)
                            {
                                if (firstIsMajorAllele)
                                {
                                    _mainmatrix.put(curind,0,((atof(tokens[i].c_str())*2) + atof(tokens[i+1].c_str())));
                                }
                                else
                                {
                                    _mainmatrix.put(curind,0,((atof(tokens[i+2].c_str())*2) + atof(tokens[i+1].c_str())));

                                }
                                for (int j = 0; j<G.phenoList.size(); j++)
                                {
                                    string xxy = G.samples[curind]._name;
                                    //     int xxyy = j;
                                    //   double xxx = G.samples[curind]._phenos[j];

                                    _mainmatrix.put(curind,j+1,G.samples[curind]._phenos[j]);
                                }
                                curind++;
                            }
                            //We have main table now - lets run through all possible combinations of phenotypes
                            //if (G.debugMode)_mainmatrix.print();
                            int _testcount = pow(2,G.phenoList.size());
                            for (int test=_testcount-1; test>=1; test--)
                            {
                                int _phenoCount = 0;
                                int _sampleCount = 0;
                                vector<bool> phenoMask = phenoMasker(test,(int)G.phenoList.size());

                                if (G.debugMode)cout << "Current mask: ";
                                for (int i = 0; i < phenoMask.size(); i++)
                                {
                                    if (phenoMask[i])
                                    {
                                        _phenoCount++;
                                        if (G.debugMode)cout<<"1";
                                    }
                                    else{if (G.debugMode)cout<<"0";}
                                }
                                if (G.debugMode)cout << endl;


                                for (int i = 0; i<_mainmatrix.getRows(); i++)
                                {
                                    bool indISOK = true;
                                    if (G.debugMode){cout << G.samples[i]._name;}
                                    if (_mainmatrix.get(i,0)==-9999){indISOK=false;}

                                    for (int j = 1; j<_mainmatrix.getCols(); j++)
                                    {

                                                //phenotypes
                                            if (phenoMask[j-1]){
                                                //if (i==0)_phenoCount++;
                                                if (_mainmatrix.get(i,j)==-9999){indISOK=false;}
                                            }
                                        if (G.debugMode){cout << " " << _mainmatrix.get(i,j);}

                                    }
                                    if (G.debugMode){cout << endl;}

                                    if (indISOK)_sampleCount++;
                                }
                                if (G.debugMode)cout << "test: " << test << " indcount: " << _sampleCount << " phenocount: " <<   _phenoCount << endl;
                                arrayD *Y = new arrayD ( _sampleCount );
                                matrixD *X = new matrixD ( _sampleCount, _phenoCount+1);
                                arrayD *W = new arrayD ( _sampleCount );
                                int curInd = 0;
                                vector <double> COPYpheno;
                                for (int i = 0; i<_mainmatrix.getRows(); i++)
                                {
                                    bool indISOK = true;
                                    if (_mainmatrix.get(i,0)==-9999){indISOK=false;}
                                    for (int j = 1; j<_mainmatrix.getCols(); j++)
                                    {

                                        if (phenoMask[j-1]){
                                            //if (i==0)_phenoCount++;
                                            if (_mainmatrix.get(i,j)==-9999){indISOK=false;}
                                        }
                                    }
                                    if (indISOK)
                                    {
                                        int z=1;

                                        //genotype
                                        Y->put(curInd, _mainmatrix.get(i,0));
                                        COPYpheno.push_back(_mainmatrix.get(i,0));

                                        W->put(curInd,1);
                                        X->put(curInd, 0,  1.0);

                                        for (int j = 1; j<_mainmatrix.getCols(); j++)
                                        {
                                            if (phenoMask[j-1])
                                            {
                                                X->put(curInd, z, _mainmatrix.get(i,j));
                                                z++;
                                            }
                                        }
                                        curInd++;
                                    }

                                }
                                if (G.debugMode)cout << "Print Y" << endl;
                                if (G.debugMode)Y->print();
                                if (G.debugMode)cout << "Print X" << endl;
                                if (G.debugMode)X->print();
                                //if (G.debugMode)
                                //{
                                //    for (int xxx=0; xxx<COPYpheno.size(); xxx++)cout<< COPYpheno[xxx];
                                //    cout <<endl;
                                //}

                                lr LR;
                                try {
                                    if (LR.lr_w(Y, X, W) && test!=_testcount)
                                    {
                                        matrixD * C = LR.covariance;
                                        C->InvertS();
                                        double x;
                                        for (int i=0; i<(unsigned int) _phenoCount;i++)
                                            for (int j=0; j<(unsigned int) _phenoCount;j++)
                                                x+=LR.Cstat->get(i)*LR.Cstat->get(j)*C->get(i,j);
                                        double testLogLikelihood = -LR.getLnLk(Y, X, W)/2;
                                        double nullLogLikelihood = -LR.nullLikelihood(COPYpheno)/2;
                                        if (G.debugMode){cout << "Linear regression done with model: " << test << endl;}
                                        for (int k = 0; k<=_phenoCount;k++)
                                        {
                                            double _pPheno = (1-studenttdistribution((_sampleCount) - _phenoCount,ddabs(LR.Cstat->get(k)/LR.SECstat->get(k))))*2;
                                            if (G.debugMode){cout << k << "\t" << (_sampleCount) << "\t" << _phenoCount << "\t" << LR.Cstat->get(k) <<  "\t" << LR.SECstat->get(k) << "\t" << _pPheno << endl;}

                                        }
                                        double likelihoodRatio = 2 * (testLogLikelihood - nullLogLikelihood);
                                        double _BIC = (-2 * testLogLikelihood) + ((_phenoCount+1) * log(_sampleCount));
                                        double _BICnull = (-2 *nullLogLikelihood) + (log(_sampleCount));
                                        //double _pModel = gsl_cdf_chisq_Q( likelihoodRatio , _phenoCount); // Overall P value - need to figure out the number of df!!!!!
                                        double _pModel;
                                        if (ddabs(likelihoodRatio)>0){_pModel = 1-chisquaredistribution(_phenoCount,ddabs(likelihoodRatio));}
                                        else {_pModel = NAN;}
                                        if (G.debugMode){cout << "Likelihood: " << testLogLikelihood << endl;}
                                        if (G.debugMode){cout << "nullLikelihood: " << nullLogLikelihood << endl;}
                                        if (G.debugMode){cout << "Model likelihood ratio: " << likelihoodRatio << endl;}
                                            if (G.debugMode){cout << "Model p: " << _pModel << endl;}
                                        if (G.debugMode){LR.covariance->print();}

                                        /*
                                        cout << chr << "\t"<<  pos << "\t" <<  markerName << "\t" << effectAllele <<"\t" << nonEffectAllele<< "\t" << _phenoCount << "\t";
                                        for (int i = 0; i < phenoMask.size(); i++) {if (phenoMask[i]){cout<<"1";}else{cout<<"0";}}
                                        cout << "\t" << LR.getLnLk(Y, X, W) <<  "\t" << -2 * nullLikelihood(COPYpheno) <<
                                        "\t" << likelihoodRatio << "\t" << _pModel << "\t" << _BIC << "\t" << _BICnull << "\t";
                                        bool phenostart = true;
                                        for (int i = 0; i < phenoMask.size(); i++) {if (phenoMask[i] && phenostart){cout<<G.phenoList[i];phenostart=false;}else if(phenoMask[i]){cout << "+" <<G.phenoList[i];}}
                                        cout << endl;
                                        */
                                        std::stringstream line, line2;
                                        if (_BIC < bestModel)
                                        {

                                            line << chr << "\t"<<  pos << "\t" <<  markerName << "\t" << effectAllele <<"\t" << nonEffectAllele<< "\t" << infoscore << "\t" << HWE(aa,aA,AA) << "\t" <<maf << "\t" << _sampleCount << "\t";
                                            if (firstIsMajorAllele){line << AA << "\t" << aA << "\t" << aa;}
                                            else {line << aa << "\t" << aA << "\t" << AA;}
                                            line << "\t" << _phenoCount << "\t";
                                            for (int i = 0; i < phenoMask.size(); i++) {if (phenoMask[i]){line<<"1";}else{line<<"0";}}
                                            line << "\t" << testLogLikelihood<<  "\t" <<  nullLogLikelihood <<
                                            "\t" << likelihoodRatio << "\t" << _pModel << "\t" << _BIC << "\t" << _BICnull << "\t";
                                            //bool phenostart = true;
                                            //for (int i = 0; i < phenoMask.size(); i++) {if (phenoMask[i] && phenostart){line<<G.phenoList[i];phenostart=false;}else if(phenoMask[i]){line << "+" <<G.phenoList[i];}}
                                            vector<string> phenosorter;
                                            for (int i = 0; i < phenoMask.size(); i++){if(phenoMask[i])phenosorter.push_back(G.phenoList[i]);}
                                            bool phenostart = true;
                                            for (int i = 0; i < phenosorter.size(); i++) {if (phenostart){line<<phenosorter[i];phenostart=false;}else{line << "+" <<phenosorter[i];}}
                                            line << "\t";
                                            sort(phenosorter.begin(),phenosorter.end());
                                            phenostart = true;
                                            for (int i = 0; i < phenosorter.size(); i++) {if (phenostart){line<<phenosorter[i];phenostart=false;}else{line << "+" <<phenosorter[i];}}

                                            if (G.printCovariance)
                                            {
                                                for (int i = 1; i < LR.Cstat->size() ;i++) line << "\t" << LR.Cstat->get(i) << "\t" <<LR.SECstat->get(i);
                                                for (int i = 1; i < LR.Cstat->size();i++)
                                                {
                                                    for (int j = i; j < LR.Cstat->size();j++)
                                                        line << "\t" << LR.covariance->get(i,j);
                                                }
                                            }

                                            line << endl;
                                            bestModelString = line.str();

                                            int k=1;
                                            for (int i = 0; i < phenoMask.size(); i++)
                                            {
                                                if (phenoMask[i])
                                                {
                                                    line2 << markerName << "\t" << effectAllele <<"\t" << nonEffectAllele<< "\t" <<  _sampleCount << "\t";
                                                    bool phenostart = true;
                                                    for (int j = 0; j < phenoMask.size(); j++)
                                                    {
                                                        if (phenoMask[j] && phenostart)
                                                        {
                                                            line2<<G.phenoList[j];
                                                            phenostart=false;
                                                        }
                                                        else if(phenoMask[j])
                                                        {
                                                            line2 << "+" <<G.phenoList[j];
                                                        }
                                                    }
                                                    line2  << "\t" << G.phenoList[i] << "\t" << LR.Cstat->get(k) << "\t" << LR.SECstat->get(k) << endl;
                                                    k++;
                                                }
                                            }
                                            bestBetasString = line2.str();
                                            bestModel = _BIC;
                                        }
                                        if (G.printAll || G.printComplex)
                                        {

                                            OUT << chr << "\t"<<  pos << "\t" <<  markerName << "\t" << effectAllele <<"\t" << nonEffectAllele<< "\t" << infoscore << "\t" << HWE(aa,aA,AA) << "\t" <<maf << "\t" << _sampleCount << "\t";
                                            if (firstIsMajorAllele){OUT << AA << "\t" << aA << "\t" << aa;}
                                            else {OUT << aa << "\t" << aA << "\t" << AA;}

                                            OUT << "\t" << _phenoCount << "\t";
                                            for (int i = 0; i < phenoMask.size(); i++) {if (phenoMask[i]){OUT<<"1";}else{OUT<<"0";}}
                                            OUT << "\t" << testLogLikelihood <<  "\t" <<  nullLogLikelihood <<
                                            "\t" << likelihoodRatio << "\t" << _pModel << "\t" << _BIC << "\t" << _BICnull << "\t";
                                            //bool phenostart = true;
                                            //for (int i = 0; i < phenoMask.size(); i++) {if (phenoMask[i] && phenostart){OUT<<G.phenoList[i];phenostart=false;}else if(phenoMask[i]){OUT << "+" <<G.phenoList[i];}}
                                            vector<string> phenosorter;
                                            for (int i = 0; i < phenoMask.size(); i++){if(phenoMask[i])phenosorter.push_back(G.phenoList[i]);}
                                            bool phenostart = true;
                                            for (int i = 0; i < phenosorter.size(); i++) {if (phenostart){OUT<<phenosorter[i];phenostart=false;}else{OUT << "+" <<phenosorter[i];}}
                                            OUT << "\t";
                                            sort(phenosorter.begin(),phenosorter.end());
                                            phenostart = true;
                                            for (int i = 0; i < phenosorter.size(); i++) {if (phenostart){OUT<<phenosorter[i];phenostart=false;}else{OUT << "+" <<phenosorter[i];}}
                                            if (G.printCovariance)
                                            {
                                                for (int i = 1; i < LR.Cstat->size() ;i++) OUT << "\t" << LR.Cstat->get(i) << "\t" <<LR.SECstat->get(i);
                                                for (int i = 1; i < LR.Cstat->size();i++)
                                                {
                                                    for (int j = i; j < LR.Cstat->size();j++)
                                                        OUT << "\t" << LR.covariance->get(i,j);
                                                }
                                            }

                                            OUT << endl;
                                            int k=1;
                                            for (int i = 0; i < phenoMask.size(); i++)
                                            {
                                                if (phenoMask[i])
                                                {
                                                    if (G.printBetas) BETAS << markerName << "\t" << effectAllele <<"\t" << nonEffectAllele<< "\t" <<  _sampleCount << "\t";
                                                    bool phenostart = true;
                                                    for (int j = 0; j < phenoMask.size(); j++)
                                                    {
                                                        if (phenoMask[j] && phenostart)
                                                        {
                                                            if (G.printBetas) BETAS<<G.phenoList[j];
                                                            phenostart=false;
                                                        }
                                                        else if(phenoMask[j])
                                                        {
                                                            if (G.printBetas) BETAS << "+" <<G.phenoList[j];
                                                        }
                                                    }
                                                    if (G.printBetas) BETAS  << "\t" << G.phenoList[i] << "\t" << LR.Cstat->get(k) << "\t" << LR.SECstat->get(k) << endl;
                                                    k++;
                                                }
                                            }
                                        }

                                        /*else if ((LR.getLnLk(Y, X, W) <= (-2 * nullLikelihood(COPYpheno))) && ((-2 * nullLikelihood(COPYpheno)) > bestModel))
                                        {
                                            line << chr << "\t"<<  pos << "\t" <<  markerName << "\t" << effectAllele <<"\t" << nonEffectAllele<< "\t" << 0 << "\t";
                                            for (int i = 0; i < phenoMask.size(); i++) {if (phenoMask[i]){line<<"0";}else{line<<"0";}}
                                            line << "\t" <<-2 * nullLikelihood(COPYpheno)  <<  "\t" << -2 * nullLikelihood(COPYpheno) <<
                                            "\t" << 1 << "\t" << 1 << "\t" << _BICnull << "\t" << _BICnull << "\tnull\n";
                                            bestModelString = line.str();
                                            bestModel = (-2 * nullLikelihood(COPYpheno));
                                        }*/
                                       // cout << setprecision(4) << test << " " << Pmodel << " "
                                       // << LR.Cstat->get(3) << " " << LR.SECstat->get(3) << " " << P1 << " "
                                       // << LR.Cstat->get(1) << " " << LR.SECstat->get(2) << " " << P2 << " "
                                       // << LR.Cstat->get(2) << " " << LR.SECstat->get(2) << " " << P3 << "\n" ;
                                    }
                                    else
                                    {
                                        if (test!=_testcount) LOG << "Collinearity problem with model: " << chr << "\t"<<  pos << "\t" <<  markerName << "\t" << test << " " << _testcount << " ";
                                        bool phenostart = true;
                                        for (int i = 0; i < phenoMask.size(); i++) {if (phenoMask[i] && phenostart){LOG<<G.phenoList[i];phenostart=false;}else if(phenoMask[i]){LOG << "+" <<G.phenoList[i];}}
                                        LOG << endl;
                                    }

                                } catch (ArgException &e)  // catch any exceptions
                                { cerr << "error: " << e.error() << " for linear model " << e.argId() << endl; }
                                delete X;
                                delete Y;
                                delete W;
                                if (G.printComplex){break;}
                            }

                            if (!G.printAll && !G.printComplex)
                            {
                                if (G.printBetas)BETAS << bestBetasString;
                                OUT << bestModelString;
                            }
                        }
                    }
                }
            }
        }
        else {cout << "Cannot read genotype file. Exit program!" << endl; exit(1);}
    }
    return true;
}
