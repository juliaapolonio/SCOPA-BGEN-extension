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

#ifndef _TOOLS_H_
#define _TOOLS_H_

#include <iostream>
#include <map>
#include <vector>
#include <stdlib.h>
#include <fstream>
#include <math.h>
#include <cctype> // std::toupper
#include <string>
#include <algorithm>
using namespace std;

void sortVec(vector <double>& x, int size);

int Tokenize(const  string& str1,
                      vector<string>& tokens,
			const string& delimiters);
string uc(string s);	//uppercase
bool checkAlleles(string & s1, string & s2);	//check if alleles are ok and change numbers to letters if necessary
string flip(string s);	//flip the alleles if 
vector <bool> phenoMasker(int, int);

string  HWE(double aa, double aA, double AA);


#endif
