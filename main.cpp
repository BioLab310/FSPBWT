
#include "FSPBWT.h"

void printHelp()
{
	std::cout << "Usage: program [options]\n";
	std::cout << "Options:\n";
	std::cout << "  -h, -H             Print this help message and exit\n";
	std::cout << "  -B, -b <value>     Set the value of B (default: 64)\n";
	std::cout << "  -F, -f <value>     Set the value of F (default: 2)\n";
	std::cout
			<< "  -i, -I <file>      Specify the panel file (default: panel.vcf)\n";
	std::cout
			<< "  -o, -O <file>      Specify the output file (default: panel_B_F_L.txt)\n";
	std::cout << "  -L, -l <value>     Set the value of L (default: 500)\n";
	std::cout
			<< "  -m, -M <mode>      in for in-panel query; out for out-panel query\n";
	std::cout << "  -q, -Q <file>      Specify the query file(default: panel.vcf)\n";
	std::cout << "  -e, -even          Set the even(position-based) mode (default: count-based)\n";
	std::cout << std::endl;
}


int main(int argc, char *argv[])
{
	int B = 64, F = 2, L = 500;
	string panelFile = "panel.vcf";
	string outputFile = "";
	string mode = "in";
	string queryFile="query.vcf";
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

