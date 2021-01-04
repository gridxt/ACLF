#include "LoadData.h"

int LoadCaseData(char *filename, double **Bus_Data, unsigned long long int *Num_of_Bus, double **Branch_Data, unsigned long long int *Num_of_Branch, double **Gen_Data, unsigned long long int *Num_of_Gen,
	double **Gen_Cost_Data, unsigned long long int *Num_of_Gen_Cost, char Bus_Data_Label[], char Branch_Data_Label[], char Gen_Data_Label[], char GenCost_Data_Label[], char BaseMVA_Label[], char Convert_Kilo_Mega_Label[], char Convert_Ohm_PU_Label[], char End_Data_Label[],
	int *BusDataColumns, int *BranchDataColumns, int *GenDataColumns, int *GenDataCostColumns, double *BaseMVA, int *Convert_Kilo, int *Convert_Ohm, double *timer)
{
	// CPU execution time to run this function //
	double start, stop;
	start = dsecnd();

	// Local variables
	FILE* fp;
	char ReadLine[bufSize]; /* Contains the line that was read from the file */
	int Bus_Data_Detected = 0; /* Set to 1 if Bus data label is detected */
	int Branch_Data_Detected = 0; /* Set to 1 if Bus data label is detected */
	int Gen_Data_Detected = 0; /* Set to 1 if Generator data label is detected */
	int GenCost_Data_Detected = 0; /* Set to 1 if Generators' Cost data label is detected */
	int Kilo_Mega_Detected = 0; /* Set to 1 if power values are in kilo and need to be converted to mega */
	int Ohm_PU_Detected = 0; /* Set to 1 if line values are in ohm and need to be converted to PU */
	int BaseMVA_Detected = 0; /* Set to 1 if BaseMVA data is detected */
	unsigned long long int i = 0, j = 0, k = 0;
	char *eptr;

	/* Open the data file */
	if ((fp = fopen(filename, "r")) == NULL)
	{
		perror("Error in opening the file!");
		return 0;
	}
//printf("\n %s", Convert_Kilo_Mega_Label);
	/* Start reading the file line-by-line */
	while (fgets(ReadLine, sizeof(ReadLine), fp) != NULL)
	{
		ReadLine[strlen(ReadLine) - 1] = '\0'; // eat the newline fgets() stores

											   /* If a known data label is detected, set the respective flag to 1 */
		if (strstr(ReadLine, Bus_Data_Label)) { Bus_Data_Detected = 1; }
		if (strstr(ReadLine, Branch_Data_Label)) { Branch_Data_Detected = 1; }
		if (strstr(ReadLine, Gen_Data_Label)) { Gen_Data_Detected = 1; }
		if (strstr(ReadLine, GenCost_Data_Label)) { GenCost_Data_Detected = 1; }
		if (strstr(ReadLine, BaseMVA_Label)) { BaseMVA_Detected = 1; }
		if (strstr(ReadLine, Convert_Kilo_Mega_Label)) { Kilo_Mega_Detected = 1; }
		if (strstr(ReadLine, Convert_Ohm_PU_Label)) { Ohm_PU_Detected = 1; }
		if (strstr(ReadLine, End_Data_Label)) { Bus_Data_Detected = 0; Branch_Data_Detected = 0; Gen_Data_Detected = 0; GenCost_Data_Detected = 0; BaseMVA_Detected = 0; Kilo_Mega_Detected = 0; Ohm_PU_Detected = 0;i = 0; j = 0; k = 0;}


 		//****** if Base MVA is detected **********
		if (BaseMVA_Detected == 1) {
			i++;
			if (i>1) { // to ignore label line and start from the next line
				char *pch;
				pch = strtok(ReadLine, " \t"); // breaks string str into a series of tokens using the delimiter delim.
				*BaseMVA = strtod(pch, NULL); // Convert the string to float
                pch = strtok(NULL, " \t"); // Remove the saved word (Bus_Data[j] from pch
				} // End of while

			}

      //****** if Kilo to Mega conversion is detected **********
		if (Kilo_Mega_Detected == 1) {
			i++;
			if (i>1) { // to ignore label line and start from the next line
				char *pch;
				pch = strtok(ReadLine, " \t"); // breaks string str into a series of tokens using the delimiter delim.
				if (strtol(pch, &eptr, 2) == 1) *Convert_Kilo = 1; // convert the string to integer and check if it is set to 1
                pch = strtok(NULL, " \t"); // Remove the saved word (Bus_Data[j] from pch
				} // End of while

			}

        //****** if Ohm to PU conversion is detected **********
		if (Ohm_PU_Detected == 1) {
			i++;
			if (i>1) { // to ignore label line and start from the next line
				char *pch;
				pch = strtok(ReadLine, " \t"); // breaks string str into a series of tokens using the delimiter delim.
				if (strtol(pch, &eptr, 2) == 1) *Convert_Ohm = 1; // convert the string to integer and check if it is set to 1
                pch = strtok(NULL, " \t"); // Remove the saved word (Bus_Data[j] from pch
				} // End of while

			}

		//****** if Buses Data is detected **********
		if (Bus_Data_Detected == 1) {
			i++;
			if (i>1) { // to ignore label line and start from the next line
				char *pch;
				pch = strtok(ReadLine, " \t"); // breaks string str into a series of tokens using the delimiter delim.

				while (pch != NULL) { // This loop save each word in pch into the Bus_Data[] array
                    k++;
					if (strstr(pch, ";")) { // Search if the last character of the string is ;
						pch[strlen(pch) - 1] = 0; // Remove ; form the word
                        *BusDataColumns = k; // number of columns in each bus data row
                        k = 0;
					}
					(*Bus_Data)[j] = strtod(pch, NULL); // Convert the string to float
					pch = strtok(NULL, " \t"); // Remove the saved word (Bus_Data[j] from pch
					j++;
					*Bus_Data = realloc(*Bus_Data, (j + 1) * sizeof(double));
					//if (j==98){Bus_Data = (double*)realloc(Bus_Data, 200 * sizeof(double));}
					*Num_of_Bus = j; // keep track of number of elements in Bus_Data[] array (not possible to do this by functions because of malloc initiation)
				} // End of while

			}

		}
		//*****************************************


		//****** if Branch Data is detected **********
		if (Branch_Data_Detected == 1) {
			i++;
			if (i>1) { // to ignore label line and start from the next line
				char *pch;
				pch = strtok(ReadLine, " \t"); // breaks string str into a series of tokens using the delimiter delim.

				while (pch != NULL) { // This loop save each word in pch into the *Branch_Data[] array
                    k++;
					if (strstr(pch, ";")) { // Search if the last character of the string is ;
						pch[strlen(pch) - 1] = 0; // Remove ; form the word
						*BranchDataColumns = k;
						k = 0;
					}
					(*Branch_Data)[j] = strtod(pch, NULL); // Convert the string to float
					pch = strtok(NULL, " \t"); // Remove the saved word (*Branch_Data[j] from pch
					j++;
					*Branch_Data = realloc(*Branch_Data, (j + 1) * sizeof(double));
					*Num_of_Branch = j; // keep track of number of elements in *Branch_Data[] array (not possible to do this by functions because of malloc initiation)
				} // End of while

			}

		}
		//*****************************************


		//****** if Generators Data is detected **********
		if (Gen_Data_Detected == 1) {
			i++;
			if (i>1) { // to ignore label line and start from the next line
				char *pch;
				pch = strtok(ReadLine, " \t"); // breaks string str into a series of tokens using the delimiter delim.

				while (pch != NULL) { // This loop save each word in pch into the Gen_Data[] array
                    k++;
					if (strstr(pch, ";")) { // Search if the last character of the string is ;
						pch[strlen(pch) - 1] = 0; // Remove ; form the word
                        *GenDataColumns = k;
                        k = 0;
					}

					(*Gen_Data)[j] = strtod(pch, NULL); // Convert the string to float
					pch = strtok(NULL, " \t"); // Remove the saved word (Gen_Data[j] from pch
					j++;
					*Gen_Data = realloc(*Gen_Data, (j + 1) * sizeof(double));
					*Num_of_Gen = j; // keep track of number of elements in Gen_Data[] array (not possible to do this by functions because of malloc initiation)
				} // End of while

			}

		}
		//*****************************************


		//****** if Generators' Cost Data is detected **********
		if (GenCost_Data_Detected == 1) {
			i++;
			if (i>1) { // to ignore label line and start from the next line
				char *pch;
				pch = strtok(ReadLine, " \t"); // breaks string str into a series of tokens using the delimiter delim.

				while (pch != NULL) { // This loop save each word in pch into the Gen_Cost_Data[] array
                    k++;
					if (strstr(pch, ";")) { // Search if the last character of the string is ;
						pch[strlen(pch) - 1] = 0; // Remove ; form the word
                        *GenDataCostColumns = k;
                        k = 0;
					}

					(*Gen_Cost_Data)[j] = strtod(pch, NULL); // Convert the string to float
					pch = strtok(NULL, " \t"); // Remove the saved word (Gen_Cost_Data[j] from pch
					j++;
					*Gen_Cost_Data = realloc(*Gen_Cost_Data, (j + 1) * sizeof(double));
					*Num_of_Gen_Cost = j; // keep track of number of elements in Gen_Cost_Data[] array (not possible to do this by functions because of malloc initiation)
				} // End of while

			}

		}
		//*****************************************



	} // End of main while loop



  fclose(fp); // close the case file

  // Calculate and return the CPU execution time
  stop = dsecnd();;
  timer[0] = stop - start;

  return 1;

} // End of LoadCasedata()
