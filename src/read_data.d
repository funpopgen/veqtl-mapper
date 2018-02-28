module read_data;

/*
This module provides functions for extracting data from bed and vcf files and storing them in appropriate structs

For phenotype data, the struct stores: the gene name, the chromosome and base pair location of the TSS and the expression values in an array

For genotype data the struct stores: the SNP ID, and the genotype values for the inidividuals, and provides a field for storing the correlation between distance and genotype for that SNP

The genotype data is extracted by calling tabix in a separate process.

There is also a function that opens a file for writing the output and writing the header.

 */
import arg_parse : Opts;
import core.stdc.stdlib : exit;
import std.algorithm : countUntil, map, max;
import std.array : array, split;
import std.conv : to;
import std.exception : enforce;
import std.format : format;
import std.process : pipeShell, Redirect, wait;
import std.range : indexed, iota;
import std.stdio : File, readln, stderr, stdout, writeln;

class InputException : Exception
{
  pure this(string s)
  {
    super(s);
  }
}

struct Phenotype
{
  string geneName;
  string chromosome;
  size_t location;
  double[] values;
  string[] eqtls;

  this(const char[] line, const size_t[] indices)
  {
    auto splitLine = line.split;
    geneName = splitLine[3].to!string;
    chromosome = splitLine[0].to!string;
    location = splitLine[2].to!size_t;
    values = splitLine[4 .. $].indexed(indices).map!(a => to!double(a)).array;
    if (countUntil!"a != b"(values, values[0]) == -1)
      throw new InputException("");
  }
}

struct Genotype
{
  string snpId;
  double[] values;
  double cor;

  this(const char[] line, const size_t[] indices, const long loc, const bool gt)
  {
    auto splitLine = line.split;
    snpId = format!"\t%-(%s\t%)\t"(splitLine[0 .. 4]);
    values = splitLine[4 .. $].indexed(indices).map!(a => getDosage(a, loc, gt)).array;
    if (countUntil!"a != b"(values, values[0]) == -1)
      throw new InputException("");
    cor = 0;
  }
}

double getDosage(const char[] field, const long loc, const bool gt)
{
  auto fieldSplit = field.split(':');
  enforce(fieldSplit.length > loc, new InputException(""));

  return gt ? cast(ubyte) fieldSplit[loc][0] + cast(ubyte) fieldSplit[loc][2] - 96
    : fieldSplit[loc].to!double;
}

auto readBed(const Opts opts)
{
  File bedFile;
  try
  {
    bedFile = File(opts.bed);
  }
  catch (Exception e)
  {
    stderr.writeln(e.msg);
    exit(1);
  }

  bedFile.readln;

  foreach (i; iota((opts.jobNumber - 1) * opts.genes))
  {
    try
    {
      bedFile.readln;
    }
    catch (Exception e)
    {
      stderr.writeln("Too few phenotypes in bed file");
      exit(1);
    }
  }

  if (bedFile.eof)
  {
    stderr.writeln("Too few phenotypes in bed file");
    exit(1);
  }

  Phenotype[] phenotype;

  auto geneCount = opts.genes == 0 ? uint.max : opts.genes;

  foreach (line; bedFile.byLine)
  {
    if (geneCount == 0)
      break;

    try
    {
      phenotype ~= Phenotype(line, opts.phenotypeLocations);
    }
    catch (InputException e)
    {
      if (opts.verbose)
      {
        stderr.writeln("Gene ", line.split[3], " is constant. Nothing to analyse.");
      }
    }
    catch (Exception e)
    {
      if (opts.verbose)
      {
        stderr.writeln("Failed to read gene ", line.split[3], ".");
      }
    }
    geneCount--;
  }

  if (phenotype.length == 0)
  {
    stderr.writeln("No phenotypes read from file.");
    exit(1);
  }

  return phenotype;
}

auto readGenotype(const Opts opts, string chrom, size_t location, string geneName)
{

  immutable size_t start = location < opts.window ? 0 : location - opts.window;

  string tabixCommand = "tabix " ~ opts.vcf ~ " " ~ chrom ~ ":" ~ start.to!string ~ "-" ~ (
      location + opts.window).to!string ~ " | grep -v '#' | cut -f1,2,4,5,10-";

  auto pipes = pipeShell(tabixCommand, Redirect.stdout);
  scope (exit)
    wait(pipes.pid);

  Genotype[] genotype;

  foreach (line; pipes.stdout.byLine)
  {
    try
    {
      genotype ~= Genotype(line, opts.genotypeLocations, opts.loc, opts.gt);
    }
    catch (Exception e)
    {
    }
  }

  if (genotype.length == 0)
  {
    stderr.writeln("Failed to extract any useful SNPs for ", geneName,
        ". Run: \n\n\"", tabixCommand, "\"\n\nfor more information.\n");
  }

  return genotype;
}

auto readCovs(const Opts opts)
{
  double[] covariates;

  File covFile;
  try
  {
    covFile = File(opts.cov);
  }
  catch (Exception e)
  {
    stderr.writeln(e.msg);
    exit(1);
  }

  covFile.readln;
  foreach (line; covFile.byLine)
  {
    covariates ~= line.split.indexed(opts.covLocations).map!(a => a.to!double).array;
  }
  return covariates;
}

auto readEqtls(string eqtlFilename, ref Phenotype[] phenotypes)
{
  string[][string] eqtlList;
  foreach (e; phenotypes)
    eqtlList[e.geneName] = [];

  File eqtlFile;
  try
  {
    eqtlFile = File(eqtlFilename);
  }
  catch (Exception e)
  {
    stderr.writeln(e.msg);
    exit(1);
  }

  foreach (line; eqtlFile.byLine)
  {
    auto splitLine = line.split;
    if (splitLine[0] in eqtlList)
    {
      eqtlList[splitLine[0].to!string] ~= format!"\t%-(%s\t%)\t"(splitLine[1 .. 5]);
    }
  }

  foreach (ref e; phenotypes)
  {
    e.eqtls = eqtlList[e.geneName].dup;
  }
}

auto makeOut(const Opts opts)
{
  File outFile;

  try
  {
    if (opts.output == "")
      outFile = stdout;
    else
      outFile = File(opts.output, "w");

    if (!opts.noheader)
      outFile.writeln("GENE\tCHROM\tPOS\tREF\tALT\tCOR\tP\tPermP\tFWER\tBETA");

  }
  catch (Exception e)
  {
    stderr.writeln("Failed to write to output file. ", e.msg);
    exit(1);
  }
  return outFile;

}
