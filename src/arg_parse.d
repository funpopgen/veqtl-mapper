module arg_parse;

import std.array : split;
import core.stdc.stdlib : exit;
import std.array : array;
import std.algorithm : canFind, countUntil, filter, map;
import std.conv : to, ConvException;
import std.exception : enforce;
import std.file : exists;
import std.process : pipeShell, Redirect, wait;
import std.range : indexed, iota;
import std.stdio : File, writeln, stderr;
import std.string : chomp;

class Opts
{
  import std.getopt;

  //write appropriate string and quit
  bool version_ = false;
  //phenotype and genotype ids are given
  bool noheader = false;
  //number of genotype columns to skip, and phenotype column
  uint genes = 0;
  uint jobNumber = 0;
  //permutation numbers and seeds
  uint[] perms;
  //file names
  string vcf = "";
  string bed = "";
  string output = "";
  size_t window = 1_000_000;
  bool gt = false;
  long loc = 0;
  bool nocheck = false;
  size_t[] genotypeLocations;
  size_t[] phenotypeLocations;

  bool het = false;

  auto processOptions(string[] args)
  {
    arraySep = ",";
    // dfmt off
    auto options = getopt(args,
			  "bed", "Phenotype file [last argument].\n", &bed,
			  "vcf", "Genotype file.\n", &vcf,
			  "out|o", "Output file [stdout].\n", &output,
			  "genes", "This specifies the number of genes to be analysed in each job.\n", &genes,
			  "job-number", "Split the analysis into a number of smaller runs which can run in parallel on a cluster. This option specifies which of the sub-analyses should be run.\n", &jobNumber,
			  "het", "Tests for variance differences in the heterozygous relative to homozygous groups, an indication of parent of origin effects.\n", &het,
			  "perm", "Specify number of permutations, with optional seed. One following number indicates the number of permutations, two comma separated numbers gives the number of permutations and the seed.\n", &perms,
			  "window", "The size in base pairs of the cis window around the transcription start site of the gene [1,000,000].\n", &window,
			  "noheader", "Suppress writing of header line.\n", &noheader,
			  "nocheck", "Do not use the header lines to match genotype and phenotype, assume samples already match.\n", &nocheck,
			  "version", "Display version information.\n", &version_,
			  );
    // dfmt on
    return options;
  }

  this(string[] args)
  {
    bool noArgs = args.length == 1;
    try
    {
      auto options = processOptions(args);
      if (options.helpWanted || noArgs)
      {
        defaultGetoptPrinter("VEQM: mapping variance associated loci.
", options.options);
        exit(0);
      }

      if (version_)
        giveHelp(versionString);

      matchIds();

      if (bed == "" && args.length > 1)
        bed = args[$ - 1];

      if (perms.length == 0)
        perms = [10_000];
    }
    catch (Exception e)
    {
      writeln(e.msg, "\n");

      auto p = ["dummy"];
      auto options = processOptions(p);
      defaultGetoptPrinter("VEQM: mapping variance associated loci.
", options.options);

      exit(1);
    }
  }

  private void matchIds()
  {

    if (!vcf.exists)
    {
      stderr.writeln("Error: genotype file ", vcf, " does not exist.");
      exit(1);
    }

    if (!(vcf ~ ".tbi").exists && !(vcf ~ ".csi").exists)
    {
      stderr.writeln("Error: Neither ", vcf, ".tbi nor ", vcf,
          ".csi files are present, meaning genotype file hasn't been indexed with tabix or bcftools.");
      exit(1);
    }

    string[] phenotypeIds;

    try
    {
      auto bedFile = File(bed);

      phenotypeIds = bedFile.readln.chomp.split[4 .. $];
    }
    catch (Exception e)
    {
      stderr.writeln("Failed to read phenotype IDs. ", e.msg);
      exit(1);
    }

    if (nocheck)
    {
      genotypeLocations = iota(phenotypeIds.length).array;
      phenotypeLocations = iota(phenotypeIds.length).array;
    }
    else
    {
      string[] genotypeIds;
      try
      {
        auto pipes = pipeShell("zcat " ~ vcf ~ " | grep -v '##' | head -2", Redirect.stdout);
        scope (exit)
          wait(pipes.pid);

        auto line = pipes.stdout.readln.chomp;

        genotypeIds = line.split[9 .. $].to!(string[]);

        auto formatField = pipes.stdout.readln.chomp.split[8].split(':');
        loc = countUntil(formatField, "DS");
        if (loc == -1)
        {
          loc = countUntil(formatField, "GT");
          gt = true;
          if (loc == -1)
          {
            stderr.writeln("DS and GT fields are both missing from the vcf file.");
            exit(1);
          }
        }
      }
      catch (Exception e)
      {
        stderr.writeln("Failed to read genotype IDs. ", e.msg);
        exit(1);
      }

      phenotypeLocations = genotypeIds.map!(a => phenotypeIds.countUntil(a))
        .filter!(a => a != -1).array.to!(size_t[]);

      genotypeLocations = iota(genotypeIds.length).filter!(
          a => phenotypeIds.canFind(genotypeIds[a])).array;

      if (phenotypeIds.indexed(phenotypeLocations)
          .array != genotypeIds.indexed(genotypeLocations).array)
      {
        stderr.writeln("Failed to match IDs. THIS SHOULD NEVER HAPPEN");
        exit(1);
      }

      if (genotypeLocations.length == 0 || phenotypeLocations.length == 0)
      {
        stderr.writeln("No individuals to analyse.");
        exit(0);
      }
    }
  }
}

static immutable string versionString = "VEQM  -  Variance association mapping for molecular phenotypes version 1.0.0";

void giveHelp(immutable string quitString)
{
  import std.compiler;

  writeln(quitString);
  writeln("Compiled with ", name, " ", version_major, ".", version_minor);
  exit(0);
}
