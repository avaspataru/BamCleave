/*************************************************************************
bamCleave.cpp 


Splits a BAM file into two bamfiles with separate headers


*************************************************************************/
#include "stringEx.h"
#include "printEx.h"
#include "libParser.h"
#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "toolkit/bamtools_sort.h"
#include<fstream>

using namespace std;
using namespace BamTools;

struct genomeStats
{
	int mappedPairs,singleMapped,chimeras;
	genomeStats():mappedPairs(0),singleMapped(0),chimeras(0){};
	void print(TsvFile & f)
	{
		f.print("Mapped pairs",mappedPairs);
		f.print("Singletons",singleMapped);
		f.print("Chimeras",chimeras);
		f.print();
	}
};



struct statistics
{
	genomeStats firstGenome,secondGenome;
	int unmapped;
	statistics() : unmapped(0){}; 

	void print(TsvFile & f)
	{
		f.print("First Genome");
		firstGenome.print(f);
		f.print("Second Genome");
		secondGenome.print(f);
		f.print();
		f.print("Unmapped reads",unmapped);
	}
	
};



class mapData : public map<string,string>
{

public:
	bool open(const string & filename)
	{
		std::ifstream f;
		stringEx line,a,b;
		f.open(filename.c_str());
		if (f.is_open())
		while (!f.eof()) {
			std::getline(f,line);
			if(line)
			{
				parseTsv(line,a,b);
				emplace(a,b);
			}
		}
		return true;
	}
	const string & mappedVal(const string & a)
	{
		const_iterator i = find(a);
		if (i != end())
			return (i->second);
		static string null;
		return null;
	}
};

struct  CellData
{
	CellData(string name) :name(name), output(NULL), count(1),saved(0) {};
	~CellData() {
		if (output) delete(output);
	};
	BamWriter * output;
	stringEx name;
	int count;
	int saved;
};

struct  bamFilesContainer : public map<string, CellData>
{
	int fileCount;
	SamHeader & header;
	RefVector & references;
	stringEx fileRoot;
	bamFilesContainer(SamHeader & header, RefVector & references, stringEx fileRoot) :header(header), references(references),
		fileRoot(fileRoot), fileCount(0) {};
	void preadd(stringEx cell = "")
	{
		iterator j = find(cell);
		if (j == end())
			emplace(cell, fileRoot);
		else
			j->second.count++;
	};


	bool add(BamAlignment & al, stringEx cell = "")
	{
		iterator j = find(cell);
		if (j->second.output)
		{
			j->second.output->SaveAlignment(al);
			j->second.saved++;
			return true;
		}
		return false;
	}
	void initialise(stringEx cell = "")
	{
		iterator j = find(cell);
		if (cell)
			j->second.name.append("_", cell, ".bam");
		else
			j->second.name.append("_all.bam");
		j->second.output = new BamWriter();
		j->second.output->Open(j->second.name, header, references);
	}
	
	void initialiseGroup(int groupID, stringEx cell = "") {
		iterator j = find(cell);
		stringEx fileName = "group_";
		fileName.append(std::to_string(groupID), ".bam");

		if (cell) {
			j->second.name = fileName;
			j->second.output = new BamWriter();
			j->second.output->Open(j->second.name, header, references);
		}
	}
	
	void closeAndIndex(BamReader & reader)
	{
		for (auto & i : This)
		{
			if (i.second.output)
			{
				i.second.output->Close();
				delete(i.second.output);
				i.second.output = NULL;
				reader.Open(i.second.name);
				reader.CreateIndex();
				reader.Close();
			}
		};
	};
};

map<string, int> getGroupTable(stringEx groupFileName) {
	//read the groupings from file and put them in a map
	map<string, int> groupTable;
	std::ifstream groupFileReader(groupFileName, std::ifstream::in);
	int cnt = 0;
	printf("Reading in the groups.\n");
	char c = groupFileReader.get();
	while (c != EOF && (c>='A' && c<='Z')) {
		string cellID = "";
		cellID += c;
		while(c!='-') { //reading the cellID
			c = groupFileReader.get();
			if(c!='-')
				cellID += c;
		}
		c = groupFileReader.get();
		int groupID = 0;
		while (c != '\n' && (c>='0'&&c<='9')) { //reading the groupID
			groupID = groupID * 10 + (c - '0');
			c = groupFileReader.get();
		}
		groupTable[cellID] = groupID;
		//printing progress
		/*if (cnt % 100 == 0) {
			printf("Have read ", cnt, " cells.\n");
			printf("Last read ", cellID, groupID, "\n");
		}
		cnt++;*/
		c = groupFileReader.get(); //skip the enter
	}
	printf("Processed cell to group mapping.\n");
	return groupTable;
}

int main(int argc, char** argv)
{
	stringEx bamFileName;
	stringEx groupFileName;
	stringEx chromosomeMapFileName;
	stringEx outputRoot;
	stringEx cellSelectFile;
	stringEx prefix;
	stringEx tagID = "XC";
	stringEx nameTag;
	int maxCells = 1000;
	bool cell = false;
	bool allReads = true;
	bool groupOption = false;

	if(argc < 1)
		exitFail("Error: parameter wrong!");
	else if(argc == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("     bamCleave                   \n");
		printf(" -b bam filename			     \n");
		printf(" -m chromosome mapping file	     \n");
		printf(" -g groups based on table file            \n");
		printf(" -p \"XYZ\"  Extracts all chromosomes beginning XYZ into a separate file or files\n");
		printf("      stripping XYZ from the chromosome names\n");
		printf(" -c <N> Creates bam files for the top N individual cells  \n");
//		printf(" -c <filename> Creates bam files for individual cells listed in <filename>  \n");
		printf(" -t <XY> Specify the tag to be used to identify single cell identity (default XC)\n");
		printf(" -n \"X\" Tag to be used to identify single cell identity in \n     read name up to character string \"X\"\n");
		printf(" -o output fileroot (default = working directory/<bam filename>)\n");

		exit(EXIT_SUCCESS);
	}

	string commandLine = argv[0];
	int ni = 1;

	while(ni < argc)
	{
		commandLine += stringEx(" ",argv[ni]);
		if(strcmp(argv[ni], "-b") == 0)
		{
			bamFileName = argv[++ni];
			commandLine += stringEx(" ", argv[ni]);
		}
		else if (strcmp(argv[ni], "-g") == 0)
		{
			groupFileName = argv[++ni];
			groupOption = true;
			commandLine += stringEx(" ", argv[ni]);
		}
		else if(strcmp(argv[ni], "-m") == 0)
		{
			chromosomeMapFileName = argv[++ni];
			allReads = false;
			commandLine += stringEx(" ", argv[ni]);
		}
		else if(strcmp(argv[ni], "-o") == 0)
		{
			outputRoot = argv[++ni];
			commandLine += stringEx(" ", argv[ni]);
		}
		else if (strcmp(argv[ni], "-p") == 0)
		{
			prefix = argv[++ni];
			allReads = false;
			commandLine += stringEx(" ", argv[ni]);
		}
		else if (strcmp(argv[ni], "-n") == 0)
		{
			nameTag = argv[++ni];
			commandLine += stringEx(" ", argv[ni]);
		}
		else if (strcmp(argv[ni], "-c") == 0)
		{
			cell = true;

			char* p;
			maxCells = strtol(argv[++ni], &p, 10);
			if (*p)
			{
				cellSelectFile = argv[ni];
				maxCells = -1;
			}
			commandLine += stringEx(" ", argv[ni]);
		}
		else if (strcmp(argv[ni], "-t") == 0)
		{
			tagID = argv[++ni];
			commandLine += stringEx(" ", argv[ni]);
		}
		ni++;
	}

#ifdef _WIN32
	if (maxCells > 500)
		_setmaxstdio(maxCells + 10);
#endif

	map<string, int> groupTable;
	set<int> groupIDs;
	if (groupOption) {
		groupTable = getGroupTable(groupFileName);
		for (std::map<string, int>::iterator it = groupTable.begin(); it != groupTable.end(); ++it) {
			if (!groupIDs.count((int)it->second))
				groupIDs.insert(it->second);
		}
		printf("Number of groups found: ", groupIDs.size(),"\n");
		printf("Number of cells read: ", groupTable.size(), "\n");
	}
	
	mapData chromosomeMap;
	if (chromosomeMapFileName)
		chromosomeMap.open(chromosomeMapFileName);


	BamReader reader;
	if (!reader.Open(bamFileName) ) 
		return EXIT_FAILURE;
 
	//	Get the header/sequence information from the existing bam file
	SamHeader header = reader.GetHeader();
	RefVector references = reader.GetReferenceData();

	//	And then create the new header/sequence information
	SamHeader header1 = header,header2 = header;
	RefVector references1, references2;

	header1.Sequences.Clear();
	header2.Sequences.Clear();

	struct destination
	{
		destination(int file, int refId) :file(file), refId(refId) {};
		int file;
		int refId;
	} ;
	struct refIdListType : public vector<destination>
	{
		bool allReads;
		refIdListType(bool allReads) : allReads(allReads) {};
		destination operator[](int i)
		{
			if (allReads)
				return destination(1, i);

			return vector<destination>::operator[](i);
		}

	} refIdList(allReads);


	int refCount = 0;
	if (allReads)
	{
		header1 = header;
		references1 = references;
		header2 = header;
		references2 = references;
	}
	else
	{
		for (SamSequenceConstIterator i = header.Sequences.Begin(); i != header.Sequences.End(); i++)
		{
			mapData::iterator chromMap = chromosomeMap.find(i->Name);
			if ((prefix) && stringEx(i->Name).startsWith(prefix))
			{
				SamSequence ss(*i);
				stringEx newName = (i->Name).substr(prefix.size());
				ss.Name = newName;
				header1.Sequences.Add(ss);
				refIdList.emplace_back(1, header1.Sequences.Size() - 1);
				references[refCount].RefName = newName;
				references1.emplace_back(references[refCount++]);
			}
			else if (chromMap != chromosomeMap.end())
			{
				SamSequence ss(*i);
				stringEx newName = chromMap->second;
				ss.Name = newName;
				header1.Sequences.Add(ss);
				refIdList.emplace_back(1, header1.Sequences.Size() - 1);
				references[refCount].RefName = newName;
				references1.emplace_back(references[refCount++]);
			}
			else
			{
				header2.Sequences.Add(*i);
				refIdList.emplace_back(2, header2.Sequences.Size() - 1);
				references2.emplace_back(references[refCount++]);
			}
		};
	};


	SamProgram thisProg;
	thisProg.CommandLine = commandLine;
	thisProg.ID = "bamCleave";

	header1.Programs.Add(thisProg);
	header2.Programs.Add(thisProg);

	BamWriter output2;

	//	Open the bam files for the human and mouse genome data
	stringEx bamFileNameCore = bamFileName.removeSuffix();


	if (outputRoot)
		bamFileNameCore = outputRoot;

	string restFileName(bamFileNameCore + "_rest.bam");

	if ((!output2.Open(restFileName, header2, references2)))
		exitFail("Failed to open second bam file");

	//	And also the log file and the file for the list of reads that were part one genome, part the other
	TsvFile chimiraFile, logFile;

	if (!allReads)
	{
		if (!chimiraFile.open(bamFileNameCore + "_chimeras.txt"))
			exitFail("Failed to open chimera file");
		chimiraFile.print("First", "", "Second", "");
	}
	if (!logFile.open(bamFileNameCore + "_split.log"))
		exitFail("Failed to open log file");


	if (prefix)
	{
		(bamFileNameCore += "_") += prefix;
		size_t lastChr = bamFileNameCore.size() - 1;
		if (strchr("-_", bamFileNameCore[lastChr]))
			bamFileNameCore = bamFileNameCore.substr(0, lastChr);
	}
	else
		bamFileNameCore += "_sel";


	bamFilesContainer bamFiles(header1, references1, bamFileNameCore);


	BamAlignment al;
	statistics statsCounts;

	int  loopCounter = 0;

	if (cell)
	{
		while (reader.GetNextAlignment(al))
		{
#ifdef _DEBUG
			if (loopCounter > 100000)
				break;
#endif
			if ((++loopCounter % 1000000) == 0)
				cout << loopCounter << " reads\n";
			destination dest = refIdList[al.RefID];
			if (dest.file == 1)
			{
				string cellId;
				if (nameTag)
				{
					cellId = al.Name.substr(0, al.Name.find_first_of(nameTag));
					bamFiles.preadd(cellId);
				}
				else
				{
					if (al.GetTag(tagID, cellId))
						bamFiles.preadd(cellId);
				}
			}
		}
		reader.Rewind();

		multimap<int, string> cellList;
		for (auto & i : bamFiles)
			cellList.emplace(i.second.count, i.first);

		printf("size of list:", groupTable.size(),"\n");
		int count = 0;
		for (multimap<int, string>::reverse_iterator i = cellList.rbegin(); (i != cellList.rend()) && (count++ <groupTable.size()); i++)
		{
			if (count % 100 == 0) {
				printf("Printed ", count, ".\n");
			}
			try {
				//bamFiles.initialise(i->second);
				bamFiles.initialiseGroup(groupTable[i->second], i->second);
			}
			catch (...)
			{
				cout << "Unable to open more than " << count << "files/n";
				bamFiles.at(i->second).output = NULL;
				break;
			}
		}

		loopCounter = 0;
	}
	else
	{
		bamFiles.preadd();
		bamFiles.initialise();
	}

	printf("Finished writing");
	if (!groupOption) {

		while (reader.GetNextAlignment(al))
		{

#ifdef _DEBUG
			if (loopCounter > 100000)
				break;
#endif

			if ((++loopCounter % 1000000) == 0)
				cout << loopCounter << " reads\n";

			destination dest = refIdList[al.RefID];

			al.RefID = dest.refId;
			//	The reference ID is within the range of chromosome identifiers associaqted with the first genome
			//	or is -1, which is an unmapped read
			if (al.MateRefID != -1)
			{
				destination dest2 = refIdList[al.MateRefID];
				if (dest.file != dest2.file)
				{
					statsCounts.firstGenome.chimeras++;
					chimiraFile.print(references[al.RefID].RefName, al.Position, chromosomeMap.mappedVal(references[al.MateRefID].RefName), al.MatePosition);
					chimiraFile.flush();
					al.SetIsMateMapped(false);
					al.MateRefID = al.RefID;
				}
				else
					al.MateRefID = dest2.refId;
			}
			if (al.IsMapped())
			{
				if (al.IsMateMapped())
					(dest.file == 1 ? statsCounts.firstGenome : statsCounts.secondGenome).mappedPairs++;
				else
					(dest.file == 1 ? statsCounts.firstGenome : statsCounts.secondGenome).singleMapped++;
			}
			else
				statsCounts.unmapped++;

			if (dest.file == 1)
			{
				if (cell)
				{
					stringEx cellId;
					if (nameTag)
					{
						cellId = al.Name.substr(0, al.Name.find_first_of(nameTag));
						if (cellId)
							if (!bamFiles.add(al, cellId))
								output2.SaveAlignment(al);
							else
								output2.SaveAlignment(al);
					}
					else
					{
						if (al.GetTag(tagID, cellId))
							if (!bamFiles.add(al, cellId))
								output2.SaveAlignment(al);
							else
								output2.SaveAlignment(al);
					}
				}
				else
					bamFiles.add(al);

			}
			else
				output2.SaveAlignment(al);
		}

		output2.Close();

		statsCounts.print(logFile);

		logFile.print("Saved cell data");
		logFile.print("Cell", "Reads", "Saved");
		for (auto & i : bamFiles)
			if (i.second.output != NULL)
				logFile.print(i.first, i.second.count, i.second.saved);
		for (auto & i : bamFiles)
			if (i.second.output == NULL)
				logFile.print(i.first, i.second.count, i.second.saved);

		cout << "Indexing first genome files" << endl;

		bamFiles.closeAndIndex(reader);

		cout << "Indexing second genome file" << endl;

		reader.Open(restFileName);
		reader.CreateIndex();
		reader.Close();
	}
}
