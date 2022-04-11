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

#ifndef __GLOBAL_H__
#define __GLOBAL_H__

#include <string>
#include <vector>
#include <map>
#include <stdlib.h>

#include "sample.h"

class global
{
public:
	global(void);
	~global(void);
	void createOutput(void);


	std::string version ;
    bool debugMode;
	int errNr;
    double threshold;
    std::vector<std::string> phenoList;                  //list of covariate column names
	std::string inputGenFile;
	std::string inputSampleFile;
    std::string inputExclFile;
	std::string outputRoot;
	std::string outputResult;
	std::string outputLog;
	std::string outputBetas;
	std::string outputError;
	std::string missingCode;
    std::map <std::string, int> exclusionList;
    bool removeMissing;
    bool printAll;
    bool printComplex;
    bool printBetas;
    bool printCovariance;
    
    int chr;
    std::vector <sample> samples;
    
};
#endif

