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

#include "global.h"

global::global(void)
{
		version = "1.0.14";
        missingCode = "NA";
        outputRoot = "scopa";
    removeMissing = false;
    printAll = false;
    printComplex = false;
    printBetas = false;
    printCovariance = false;
        threshold=0.95;
    chr=0;
}

global::~global(void)
{
}

void
global::createOutput()
{
	outputResult = outputRoot + ".result";
	outputLog = outputRoot + ".log";
    outputBetas = outputRoot + ".betas";
	outputError = outputRoot + ".err";
}

// 1.0.9 infoscore fix
// 1.0.10 using threshold for minimum info score now
// 1.0.11 some gzipped files have odd number of columns - using n-1 now to make sure that "\n" doesnt make it to read phantom sample
// 1.0.12 fixed allele flip issue - effect allele is now minor allele, which is necessary for meta-analysis. also allowed longer allele names
// 1.0.13 renamed software
