# Getting PyExoCross Program

## Download PyExoCross

Download PyExoCross from [GitHub](https://github.com/nichollsh/PyExoCross.git "GitHub").

```bash
git clone https://github.com/nichollsh/PyExoCross.git
```

## Install Python packages

```bash
pip install -Ue .
```

## Run PyExoCross

Prepare an input file *filename.inp* (see examples in the 'input' folder on GitHub) and run the program with command:

```bash
python run.py -p input_filepath
```

*Example*

If the input filepath is `/home/username/PyExoCross/input/H2O_exomol.inp`.

```bash
python run.py -p ./input/H2O_exomol.inp
# OR 
python run.py -p /home/username/PyExoCross/input/H2O_exomol.inp
```

If you want to run program in conda environment which is named as 'exomol', please use command:

```bash
/home/username/anaconda3/envs/exomol/bin/python run.py -p ./input/H2O_exomol.inp
```

If you need to run program in background, please use command:

```bash
# Discard redirected terminal output. LogFilePath still controls log-file output.
nohup python -u run.py -p ./input/H2O_exomol.inp > /dev/null 2>&1
# Save redirected terminal output separately
nohup python -u run.py -p ./input/H2O_exomol.inp > ./output/H2O_exomol.out 2>&1 &
# OR 
nohup /home/username/anaconda3/envs/exomol/bin/python -u run.py -p ./input/H2O_exomol.inp > ./output/H2O_exomol.out 2>&1 &
```

To disable PyExoCross log-file output, set this in the input file:

```text
LogFilePath                             None
```

This is independent of terminal output. Use `Verbose False` to hide normal
terminal output and progress bars.

## Resuming an interrupted run

If a cross-section run is interrupted (killed, timed out), it can be resumed with `-r` /
`--resume`. This skips (T, P) cross-section grid points whose `.xsec` output file already
exists in the output folder, except the 2 most recently modified existing files, which are
always redone (they may have been half-written when the previous run was interrupted). With
`-r` / `--resume` omitted (default), behaviour is unchanged: the output folder is cleared and
everything is recomputed.

```bash
python run.py -p ./input/H2O_exomol.inp -r
# OR
python run.py -p ./input/H2O_exomol.inp --resume
```

The `pyexocross` console-script entry point supports the same flag:

```bash
pyexocross -p ./input/H2O_exomol.inp --resume
```

`resume` has no `.inp` keyword equivalent -- it is CLI-only, or settable as a Python-API kwarg
(`resume=True`) for `px.cross_sections()`-style calls. See
[inp_mapping.md](python_api/inp_mapping.md) for details.

## Notes for input file

All information can be written in the input file. Just change the information you will use.You don't need to change any other unnecessary information.Please do not change the first column strings.
