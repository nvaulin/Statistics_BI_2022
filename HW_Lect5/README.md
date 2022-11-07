# DGE analyser

> *by* Nikita Vaulin, Skoltech <br />
> Nikita.Vaulin@skoltech.ru

This script allows you to analyze the differential genes expressions recorded in the form of tabular data.

The perfomance was tested with Python 3.9.10

## Dependencies

To use the script, install the following dependencies.

```commandline
pip install -r requirements.txt
```

## Installation and usage

To install `DGE_analysis.py` script simply download it.

For Linux:

```bash
wget -c https://raw.githubusercontent.com/nvaulin/Statistics_BI_2022/HW_L5/HW_Lect5/DGE_analysis.py
```

This script has several input arguments:

- _Input_1_ and _Input_2_ : two input .csv tables you want compare with differend genes as columns and their expression
  values
- _output_ : prefix for output .csv file name

To run the script paste in your terminal:

```commandline
python ./DGE_analysis.py <input_1> <input_2> <output>
```

It will provide you with output data containing statistical comparison of to data sets and, moreover, it will ask you
weather you like to get some plots on the certain genes to illustrate the analysis.

## Uninstallation

To uninstall this file, simply delete it as a regular file from your computer. For example, for Linux:

```bash
rm DGE_analysis.py
```

### If any questions

Do not hesitate to write me by email Nikita.Vaulin@skoltech.ru
