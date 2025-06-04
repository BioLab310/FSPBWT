//============================================================================
// Name        : FSPBWT_variable.cpp
// Author      :
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include "FSPBWT.h"

void printHelp()
{
	std::cout << "Usage: program [options]\n";
	std::cout << "Options:\n";
	std::cout << "  -h, -H             Print this help message and exit\n";
	std::cout << "  -B, -b <value>     Set the value of B (default: 128)\n";
	std::cout << "  -F, -f <value>     Set the value of F (default: 1)\n";
	std::cout
			<< "  -i, -I <file>      Specify the panel file (default: panel.vcf)\n";
	std::cout
			<< "  -o, -O <file>      Specify the output file (default: panel_B_F_L.txt)\n";
	std::cout << "  -L, -l <value>     Set the value of L (default: 20000)\n";
	std::cout
			<< "  -m, -M <mode>      in for in-panel query; out for out-panel query\n";
	std::cout << "  -q, -Q <file>      Specify the query file\n";
	std::cout
			<< "  -save, -SAVE <file>   Save the panel to the specified file\n";
	std::cout
			<< "  -load, -LOAD <file>   Load the panel from the specified file\n";
	// 添加了 -e 和 -even 选项的说明
	std::cout << "  -e, -even          Set the even mode (default: false)\n";
	std::cout << std::endl;
}

int savePanel(int B, int F, string save, string panel, bool even)
{
	if (save == "")
	{
		save += panel + "_saveFile_" + std::to_string(B) + "_"
				+ std::to_string(F) + ".bin";
	}
	if (B == 64)
	{
		FSPBWT<unsigned long long> CRY;
		CRY.F=F;
		CRY.B = B;
		CRY.T = pow(2, F); // 计算T的值，即2的F次方
		CRY.minSiteL = B * 2 - 1;
		// 构造函数用于初始化B、F、T和minSiteL的值

		int a = CRY.readVCF(panel);
		std::cout << "read panel file done: " << a << endl;
		if (even == true)
		{
			int b = CRY.makeFuzzyPanelEvenly();
			std::cout << "make fuzzy panel evenly done: " << b << endl;
		}
		else
		{
			int b = CRY.makeFuzzyPanelGlobally();
			std::cout << "make fuzzy panel globally done: " << b << endl;
		}

		int c = CRY.save(save);
		std::cout << "save panel done: " << c << endl;
		return 0;
	}
	else if (B == 128)
	{
		FSPBWT<__uint128_t> CRY;
		CRY.F=F;
		CRY.B = B;
		CRY.T = pow(2, F); // 计算T的值，即2的F次方
		CRY.minSiteL = B * 2 - 1;
		// 构造函数用于初始化B、F、T和minSiteL的值

		int a = CRY.readVCF(panel);
		std::cout << "read panel file done: " << a << endl;
		if (even == true)
		{
			int b = CRY.makeFuzzyPanelEvenly();
			std::cout << "make fuzzy panel evenly done: " << b << endl;
		}
		else
		{
			int b = CRY.makeFuzzyPanelGlobally();
			std::cout << "make fuzzy panel globally done: " << b << endl;
		}

		int c = CRY.save(save);
		std::cout << "save panel done: " << c << endl;
		return 0;
	}
	return 1;
}

int loadPanel(int B, int F, string save, string query, string output, int L,
		string mode, bool even)
{

	if (B == 64)
	{
		FSPBWT<unsigned long long> CRY;
		CRY.F=F;
		CRY.B = B;
		CRY.T = pow(2, F); // 计算T的值，即2的F次方
		CRY.minSiteL = B * 2 - 1;
		// 构造函数用于初始化B、F、T和minSiteL的值
		int a = CRY.load(save.c_str());
		cout << "load panel done " << a << endl;
		cout << save.c_str() << endl;
		if (mode == "in")
		{
			if (output == "")
			{
				output += save + "_inPanelQuery_" + std::to_string(B) + "_"
						+ std::to_string(F) + "_" + std::to_string(L) + ".txt";
			}
			int b = CRY.inPanelLongMatchQuery(L, output);
			cout << " in-panel query done " << b << endl;
		}
		else if (mode == "out")
		{
			if (output == "")
			{
				output += save + "_outPanelQuery_" + std::to_string(B) + "_"
						+ std::to_string(F) + "_" + std::to_string(L) + ".txt";
			}
			int c = CRY.readQueryVCF(query);
			cout << "read query done " << c << endl;
			int d = CRY.outPanelLongMatchQuery(L, output, even);
			cout << " out-panel query done " << d << endl;
		}
	}
	else if (B == 128)
	{
		FSPBWT<__uint128_t> CRY;
		CRY.F=F;
		CRY.B = B;
		CRY.T = pow(2, F); // 计算T的值，即2的F次方
		CRY.minSiteL = B * 2 - 1;
		// 构造函数用于初始化B、F、T和minSiteL的值
		int a = CRY.load(save.c_str());
		cout << "load panel done " << a << endl;
		if (mode == "in")
		{
			int b = CRY.inPanelLongMatchQuery(L, output);
			cout << " in-panel query done " << b << endl;
		}
		else if (mode == "out")
		{
			int c = CRY.readQueryVCF(query);
			cout << "read query done " << c << endl;
			int d = CRY.outPanelLongMatchQuery(L, output, even);
			cout << " out-panel query done " << d << endl;
		}
	}
	return 1;
}

int inPanel(int B, int F, string panel, string output, int L, bool even)
{
	if (B == 64)
	{
		FSPBWT<unsigned long long> CRY;
		CRY.F=F;
		CRY.B = B;
		CRY.T = pow(2, F); // 计算T的值，即2的F次方
		CRY.minSiteL = B * 2 - 1;
		// 构造函数用于初始化B、F、T和minSiteL的值

		int a = CRY.readVCF(panel);
		std::cout << "read panel file done: " << a << endl;
		if (even == true)
		{
			int b = CRY.makeFuzzyPanelEvenly();
			std::cout << "make fuzzy panel evenly done: " << b << endl;
		}
		else if(even == false)
		{
			int b = CRY.makeFuzzyPanelGlobally();
			std::cout << "make fuzzy panel globally done: " << b << endl;

		}

		if (output == "")
		{
			if(even)
			{
				output += "evenFuzzy_";
			}
			output += panel + "_inPanelQuery_"  + std::to_string(B)
					+ "_" + std::to_string(F) + "_" + std::to_string(L)
					+ ".txt";
		}
		int c = CRY.inPanelLongMatchQuery(L, output);
		std::cout << "in-panel query done: " << c << endl;

		string informationFile ="information_";
		if(even)
		{
			informationFile += "evenFuzzy_";
		}
		informationFile += panel + "_inPanelQuery_" + std::to_string(B) + "_"
				+ std::to_string(F) + "_" + std::to_string(L) + ".txt";
		CRY.outputInformationToFile(informationFile, "in");
		return 0;
	}
	else if (B == 128)
	{
		FSPBWT<__uint128_t> CRY;
		CRY.F=F;
		CRY.B = B;
		CRY.T = pow(2, F); // 计算T的值，即2的F次方
		CRY.minSiteL = B * 2 - 1;
		// 构造函数用于初始化B、F、T和minSiteL的值

		int a = CRY.readVCF(panel);
		std::cout << "read panel file done: " << a << endl;
		if (even == true)
		{
			int b = CRY.makeFuzzyPanelEvenly();
			std::cout << "make fuzzy panel evenly done: " << b << endl;
		}
		else if(even ==false)
		{
			int b = CRY.makeFuzzyPanelGlobally();
			std::cout << "make fuzzy panel globally done: " << b << endl;

		}
		if (output == "")
		{
			if(even)
			{
				output += "evenFuzzy_";
			}
			output += panel + "_inPanelQuery_"  + std::to_string(B)
					+ "_" + std::to_string(F) + "_" + std::to_string(L)
					+ ".txt";
		}
		int c = CRY.inPanelLongMatchQuery(L, output);
		std::cout << "in-panel query done: " << c << endl;

		string informationFile = "information_";
		if(even)
		{
			informationFile += "evenFuzzy_";
		}
		informationFile += panel + "_inPanelQuery_" + std::to_string(B) + "_"
				+ std::to_string(F) + "_" + std::to_string(L) + ".txt";
		CRY.outputInformationToFile(informationFile, "in");
		return 0;
	}
	return 1;
}

int outPanel(int B, int F, int L, string panel, string query, string output,
		bool even)
{
	if (B == 64)
	{
		FSPBWT<unsigned long long> CRY;
		CRY.F=F;
		CRY.B = B;
		CRY.T = pow(2, F); // 计算T的值，即2的F次方
		CRY.minSiteL = B * 2 - 1;
		// 构造函数用于初始化B、F、T和minSiteL的值
			int a = CRY.readVCF(panel);

		std::cout << "read panel file done: " << a << endl;
		if (even == true)
		{
			int b = CRY.makeFuzzyPanelEvenly();
			std::cout << "make fuzzy panel evenly done: " << b << endl;
		}
		else
		{
			int b = CRY.makeFuzzyPanelGlobally();
			std::cout << "make fuzzy panel globally done: " << b << endl;
		}
		int c = CRY.readQueryVCF(query);
		std::cout << "read query file done: " << c << endl;


		if (output == "")
		{
			if(even)
			{
				output += "evenFuzzy_";
			}
			output += panel + "_outPanelQuery_"  + std::to_string(B)
					+ "_" + std::to_string(F) + "_" + std::to_string(L)
					+ ".txt";
		}
		int d = CRY.outPanelLongMatchQuery(L, output, even);
		cout << "out-panel query done: " << d << endl;
		string informationFile = "information_";
		if(even)
		{
			informationFile+="enven_";
		}
		informationFile += panel + "_outPanelQuery_" + std::to_string(B) + "_"
				+ std::to_string(F) + "_" + std::to_string(L) + ".txt";
		CRY.outputInformationToFile(informationFile, "out");

		return 0;
	}
	else if (B == 128)
	{
		FSPBWT<__uint128_t> CRY;
		CRY.F=F;
		CRY.B = B;
		CRY.T = pow(2, F); // 计算T的值，即2的F次方
		CRY.minSiteL = B * 2 - 1;
		// 构造函数用于初始化B、F、T和minSiteL的值
		int a = CRY.readVCF(panel);
		std::cout << "read panel file done: " << a << endl;
		if (even == true)
		{
			int b = CRY.makeFuzzyPanelEvenly();
			std::cout << "make fuzzy panel evenly done: " << b << endl;
		}
		else
		{
			int b = CRY.makeFuzzyPanelGlobally();
			std::cout << "make fuzzy panel globally done: " << b << endl;
		}

		int c = CRY.readQueryVCF(query);
		std::cout << "read query file done: " << c << endl;

		if (output == "")
		{
			if(even)
			{
				output += "evenFuzzy_";
			}
			output += panel + "_outPanelQuery_"  + std::to_string(B)
					+ "_" + std::to_string(F) + "_" + std::to_string(L)
					+ ".txt";
		}
		int d = CRY.outPanelLongMatchQuery(L, output, even);
		cout << "out-panel query done: " << d << endl;

		string informationFile = "information_";
		if(even)
		{
			informationFile+="even_";
		}
		informationFile += panel + "_outPanelQuery_" + std::to_string(B) + "_"
				+ std::to_string(F) + "_" + std::to_string(L) + ".txt";
		CRY.outputInformationToFile(informationFile, "out");
		return 0;
	}
}

int main(int argc, char *argv[])
{
	int B = 64, F = 1, L = 130;
	string panelFile = "chr16.vcf";
	string outputFile = "";
	string mode = "in";
	string queryFile="chr22_query.vcf";
	bool save = false;
	string saveFile = "";
	bool load = false;
	string loadFile = "";
	bool even = false;

	if (argc == 1)
	{
		printHelp();
		return 0;
	}

	for (int i = 1; i < argc; i++)
	{
		string arg = argv[i];
		if (arg == "-h" || arg == "-H")
		{
			printHelp(); // 打印帮助信息
			exit(0); // 退出程序
		}
		if (arg == "-B" || arg == "-b")
		{
			if (i + 1 < argc)
			{
				B = atoi(argv[i + 1]);
				i++;  // 跳过下一个参数
			}
		}
		else if (arg == "-F" || arg == "-f")
		{
			if (i + 1 < argc)
			{
				F = atoi(argv[i + 1]);
				i++;  // 跳过下一个参数
			}
		}
		else if (arg == "-i" || arg == "-I")
		{
			if (i + 1 < argc)
			{
				panelFile = argv[i + 1];
				i++;  // 跳过下一个参数
			}
		}
		else if (arg == "-o" || arg == "-O")
		{
			if (i + 1 < argc)
			{
				outputFile = argv[i + 1];
				i++;  // 跳过下一个参数
			}
		}
		else if (arg == "-L" || arg == "-l")
		{
			if (i + 1 < argc)
			{
				L = atoi(argv[i + 1]);
				i++;  // 跳过下一个参数
			}
		}
		else if (arg == "-m" || arg == "-M")
		{
			if (i + 1 < argc)
			{
				mode = argv[i + 1];
				i++;  // 跳过下一个参数
			}
		}
		else if (arg == "-q" || arg == "-Q")
		{
			if (i + 1 < argc)
			{
				queryFile = argv[i + 1];
				i++;  // 跳过下一个参数
			}
		}

		else if (arg == "-save" || arg == "-SAVE")
		{
			save = true;
			if (i + 1 < argc)
			{
				saveFile = argv[i + 1];
				i++;  // 跳过下一个参数
			}
		}
		else if (arg == "-load" || arg == "-LOAD")
		{
			load = true;
			if (i + 1 < argc)
			{
				loadFile = argv[i + 1];
				i++;  // 跳过下一个参数
			}
		}
		else if (arg == "-e" || arg == "-even")
		{
			even = true; // 如果出现 -e 或 -even 选项，则将 even 设置为 true
		}
	}




	if (F<=0 || F>4) {
		cout << "wrong fuzzy size! must be 1/2/3/4" << endl;
	}
	if (outputFile=="")
	{
		if (mode == "out") {
			outputFile = "FSPBWT_outPanel_" + panelFile + "_" + queryFile+ "_" + std::to_string(B)
					+ "_" + std::to_string(F) + "_" + std::to_string(L)
					+ ".txt";
		}
		else if (mode == "in") {
			outputFile = "FSPBWT_inPanel_" + panelFile + "_" + std::to_string(B)
					+ "_" + std::to_string(F) + "_" + std::to_string(L)
					+ ".txt";
		}

	}
	string informationFile="Inf_FSPBWT_";
	if (mode == "out") {
		informationFile += panelFile + "_outPanelQuery_" + std::to_string(B) + "_"
				+ std::to_string(F) + "_" + std::to_string(L) + ".txt";
	}
	else if (mode == "in") {
		informationFile += panelFile + "_inPanelQuery_" + std::to_string(B) + "_"
					+ std::to_string(F) + "_" + std::to_string(L) + ".txt";

	}
	else {
		std::cout << "wrong mode! must be in / out" << endl;
	}

	// 输出参数
	cout << "Parameters:" << endl;
	cout << "  B: " << B << endl;
	cout << "  F: " << F << endl;
	cout << "  L: " << L << endl;
	cout << "  Input file: " << panelFile << endl;
	cout << "  Output file: " << outputFile << endl;
	cout << "  Mode: " << mode << endl;
	cout << "  Query file: " << queryFile << endl;
	cout << "  Save: " << (save ? "true" : "false") << endl;
	cout << "  Save file: " << saveFile << endl;
	cout << "  Load: " << (load ? "true" : "false") << endl;
	cout << "  Load file: " << loadFile << endl;
	cout << "  Even: " << (even ? "true" : "false") << endl;

	if (B==64) {
		FSPBWT<unsigned long long> CRY;
		CRY.F=F;
		CRY.B = B;
		CRY.T = pow(2, F); // 计算T的值，即2的F次方
		CRY.minSiteL = B * 2 - 1;
		// 构造函数用于初始化B、F、T和minSiteL的值
		int a = CRY.readVCF(panelFile);
		// int a = CRY.readTXT("sites.txt");
			std::cout << "read panel file done: " << a << endl;
			//1846144
		int b;
		if (even==true)
		{
			b = CRY.makeFuzzyPanelEvenly();
		}
		else
		{
			b = CRY.makeFuzzyPanelGlobally();

		}
			//2139264
			std::cout << "make fuzzy panel  done: " << b << endl;
			if (mode == "out") {
				int c = CRY.readQueryVCF(queryFile);
				cout << "read query done: " << c << endl;

				// 2146176
				int d = CRY.outPanelLongMatchQuery(L, outputFile, even);
				cout << "out-panel query done: " << d << endl;
				// 2146688
				CRY.outputInformationToFile(informationFile, "out");
			}
			else if (mode == "in") {
				int c = CRY.inPanelLongMatchQuery(L, outputFile);
				std::cout << "in-panel query done: " << c << endl;
				CRY.outputInformationToFile(informationFile, "in");
			}
			return 0;
	}
	else if (B==128) {
		FSPBWT<unsigned __uint128_t> CRY;
		CRY.F=F;
		CRY.B = B;
		CRY.T = pow(2, F); // 计算T的值，即2的F次方
		CRY.minSiteL = B * 2 - 1;
		// 构造函数用于初始化B、F、T和minSiteL的值
		int a = CRY.readVCF(panelFile);

			std::cout << "read panel file done: " << a << endl;
		int b;
		if (even==true)
		{
			b = CRY.makeFuzzyPanelEvenly();
		}
		else
		{
			b = CRY.makeFuzzyPanelGlobally();

		}
			std::cout << "make fuzzy panel  done: " << b << endl;

			if (mode == "out") {
				int c = CRY.readQueryVCF(queryFile);
				cout << "read query done: " << c << endl;
				int d = CRY.outPanelLongMatchQuery(L, outputFile, even);
				cout << "out-panel query done: " << d << endl;
				CRY.outputInformationToFile(informationFile, "out");
			}
			else if (mode == "in") {
				int c = CRY.inPanelLongMatchQuery(L, outputFile);
				std::cout << "in-panel query done: " << c << endl;
				CRY.outputInformationToFile(informationFile, "in");
			}
			return 0;

	}
	else {
		std::cout << "wrong Syllable size! s must be 64 / 128" << endl;
	}
	return 0;
}
