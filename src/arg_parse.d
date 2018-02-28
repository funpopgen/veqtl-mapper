module arg_parse;

/*

This file contains the routines for parsing the command line arguments.

It also performs the following task:

1) Checks if tabix is installed and fails nicely if not.
2) Reads the phenotype and genotype IDs and calculates the permutations required to match them.
3) Checks if the DS or GT field is present in the vcf file, using them prefereably in that order.

 */
import core.stdc.stdlib : exit;
import std.algorithm : canFind, countUntil, filter, joiner, map, setDifference,
  sort;
import std.array : array, split;
import std.conv : to;
import std.file : exists;
import std.getopt : arraySep, defaultGetoptPrinter, getopt;
import std.process : executeShell, pipeShell, Redirect, wait;
import std.range : indexed, iota;
import std.stdio : File, stderr, writeln;
import std.string : chomp;

class Opts
{

  //write appropriate string and quit
  bool version_ = false;
  bool verbose = false;
  //phenotype and genotype ids are given
  bool noheader = false;
  //number of genotype columns to skip, and phenotype column
  uint genes = 0;
  uint jobNumber = 0;
  //permutation numbers and seeds
  uint[] perms;
  bool normal = false;
  //file names
  string vcf = "";
  string bed = "";
  string output = "";
  string eqtl = "";
  string cov = "";
  size_t window = 1_000_000;
  bool gt = false;
  long loc = 0;
  bool nocheck = false;
  size_t[] genotypeLocations;
  size_t[] phenotypeLocations;
  size_t[] covLocations;

  bool het = false;

  auto processOptions(string[] args)
  {
    arraySep = ",";
    // dfmt off
    auto options = getopt(args,
			  "bed", "Phenotype file [last argument].\n", &bed,
			  "vcf", "Genotype file.\n", &vcf,
			  "out|o", "Output file [stdout].\n", &output,
			  "verbose", "Print additional information.\n", &verbose,
			  "genes", "This specifies the number of genes to be analysed in each job.\n", &genes,
			  "job-number", "Split the analysis into a number of smaller runs which can run in parallel on a cluster. This option specifies which of the sub-analyses should be run.\n", &jobNumber,
			  "het", "Tests for variance differences in the heterozygous relative to homozygous groups, an indication of parent of origin effects.\n", &het,
			  "perm", "Specify number of permutations, with optional seed. One following number indicates the number of permutations, two comma separated numbers gives the number of permutations and the seed.\n", &perms,
			  "normal", "Map the phenotype onto a normal distribution before analysis.\n", &normal,
			  "eqtl", "Control for a set of mapped eQTLs to remove v-eQTL caused by haplotype effects.\n", &eqtl,
			  "cov", "Control for covariates.\n", &cov,
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
    immutable bool noArgs = args.length == 1;
    try
    {
      auto options = processOptions(args);
      if (options.helpWanted || noArgs)
      {
        defaultGetoptPrinter("veqtl-mapper: mapping variance associated loci.
", options.options);
        exit(0);
      }

      if (version_)
        giveHelp(versionString);

      immutable auto checkTabix = executeShell("command -v tabix");

      if (checkTabix.status != 0)
      {
        stderr.writeln("Error: tabix is not installed.");
        exit(1);
      }

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
      defaultGetoptPrinter("veqtl-mapper: mapping variance associated loci.
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
      if (verbose)
      {
        stderr.writeln(phenotypeIds.length, " individuals present in phenotype file.");
      }
    }
    catch (Exception e)
    {
      stderr.writeln("Failed to read phenotype IDs. ", e.msg);
      exit(1);
    }

    if (nocheck)
    {
      if (verbose)
      {
        stderr.writeln("Assuming same individuals in genotype file.");
      }
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
        if (verbose)
        {
          stderr.writeln(genotypeIds.length, " individuals present in genotype file.");
        }

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

      if (genotypeLocations.length == 0 || phenotypeLocations.length == 0)
      {
        stderr.writeln("No individuals to analyse.");
        exit(1);
      }

      if (phenotypeIds.indexed(phenotypeLocations)
          .array != genotypeIds.indexed(genotypeLocations).array)
      {
        stderr.writeln("Failed to match IDs. THIS SHOULD NEVER HAPPEN");
        exit(1);
      }

      if (verbose && genotypeLocations.length != genotypeIds.length)
      {
        stderr.writeln(genotypeIds.indexed(setDifference(iota(genotypeIds.length),
            genotypeLocations)).joiner(", "), " dropped from genotype file.");
      }

      if (verbose && phenotypeLocations.length != phenotypeIds.length)
      {
        stderr.writeln(phenotypeIds.indexed(setDifference(iota(phenotypeIds.length),
            phenotypeLocations.dup.sort!())).joiner(", "), " dropped from phenotype file.");
      }

      if (cov != "")
      {
        string[] covIds;

        try
        {
          covIds = File(cov).readln.chomp.split;
          if (verbose)
          {
            stderr.writeln(covIds.length, " individuals present in covariates file.");
          }

          auto temp = genotypeIds.indexed(genotypeLocations).map!(a => covIds.countUntil(a)).array;
          if (temp.canFind(-1))
          {
            stderr.writefln("Individuals in genotype and phenotype file, but not in covariate file. Missing individuals are %-(%s, %).",
                genotypeIds.indexed(genotypeLocations).filter!(a => !covIds.canFind(a)));
            exit(1);
          }

          covLocations = temp.to!(size_t[]);

          if (verbose && covLocations.length != covIds.length)
          {
            stderr.writefln("%-(%s, %) dropped from covariates file.",
                covIds.indexed(setDifference(iota(covIds.length), covLocations.dup.sort!())));
          }
        }
        catch (Exception e)
        {
          stderr.writeln("Failed to read covariate IDs. ", e.msg);
          exit(1);
        }
      }
    }
  }
}

static immutable string versionString = "veqtl-mapper  -  Variance association mapping for molecular phenotypes version 1.0.2";
static immutable string commitString = chomp(cast(string) import("commit"));

void giveHelp(immutable string quitString)
{
  import std.compiler : name, version_major, version_minor;

  static string[] dateString = __DATE__.split;
  writeln(quitString, "-", commitString);
  writeln("Compiled with ", name, " ", version_major, ".", version_minor,
      " at ", __TIME__, ", ", dateString[1], " ", dateString[0], " ", dateString[2], ".");
  exit(0);
}
