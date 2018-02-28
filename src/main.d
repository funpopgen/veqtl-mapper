/* The GPL v3 License

   Copyright (C) 2016 University of Geneva.
   #
   # Author: Andrew Brown <andrew.brown@unige.ch>

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

/*

This module wraps the other modules, calling the functions:

1. parse the command line arguments
2. Read the phenotype data for the correct job chunk
3. Open a file to write the results to.
4. Call the routines to run the analysis.

It also contains the unit tests for runnning 2 sample analyses (standard
and parent of origin).

*/

import arg_parse : Opts;
import calculation : genPerms;
import core.stdc.stdlib : exit;
import read_data : makeOut, readBed, readCovs, readEqtls;
import run_analysis : analyseData;
import std.conv : to;
import std.range : enumerate;
import std.stdio : File, stderr, writeln;

version (STATICLINKED)
{
  pragma(msg, "Statically linked");
}
else
{
  pragma(lib, "gsl");
  pragma(lib, "gslcblas");
}

version (unittest)
  void main()
{
  writeln("All unit tests completed successfully.");
}

else
  void main(string[] args)
{
  pragma(msg, "veqtl-mapper");

  const auto opts = new Opts(args.to!(string[]));

  if (opts.verbose)
  {
    stderr.writeln("Looking for ", opts.het ? "parent of origin" : "variance", " effects.");
  }

  auto phenotype = readBed(opts);

  if (opts.eqtl != "")
  {
    readEqtls(opts.eqtl, phenotype);
  }

  if (opts.verbose)
  {
    stderr.writeln("Read ", phenotype.length, " phenotypes.");
  }
  auto permutations = genPerms(opts, phenotype[0].values.length);

  auto outFile = makeOut(opts);

  auto orderBuffer = new size_t[](phenotype[0].values.length);

  double[] covariates;

  if (opts.cov != "")
  {
    covariates = readCovs(opts);
  }

  foreach (ref e; phenotype.enumerate)
  {
    if (opts.verbose)
    {
      stderr.writeln("Analysing gene ", e[1].geneName, " (", e[0] + 1,
          " out of ", phenotype.length, "). TSS is ", e[1].chromosome, ":", e[1].location, ".");
    }
    analyseData(e[1], permutations, outFile, opts, orderBuffer, covariates);
  }
}

@system unittest
{
  import core.stdc.stdlib : exit;
  import std.array : split;
  import std.digest.sha : SHA1, toHexString;
  import std.file : exists, remove;
  import std.range : put;
  import std.stdio : stderr;
  import std.uuid : randomUUID;

  auto testFile = randomUUID.toString;

  while (testFile.exists)
  {
    testFile = randomUUID.toString;
  }

  scope (exit)
  {
    if (testFile.exists)
      testFile.remove;
  }

  // ./bin/veqtl-mapper --bed data/phenotype.bed --job-number 1 --genes 10 --vcf data/genotype.vcf.gz --perm 10000,4

  const auto opts = new Opts(("./bin/veqtl-mapper --bed data/phenotype.bed --job-number 1 --genes 10 --vcf data/genotype.vcf.gz --perm 10000,4 --out " ~ testFile)
      .split);

  auto phenotype = readBed(opts);

  auto permutations = genPerms(opts, phenotype[0].values.length);

  auto outFile = makeOut(opts);

  auto orderBuffer = new size_t[](phenotype[0].values.length);

  foreach (ref e; phenotype)
  {
    analyseData(e, permutations, outFile, opts, orderBuffer);
  }

  outFile.close;

  SHA1 hash;
  hash.start;
  put(hash, File(testFile).byChunk(1024));
  assert(toHexString(hash.finish) == "26F6F2F43BB3543FE7A8E060F0F9A088A74EACA1");
  stderr.writeln("Passed: variance test.");

  // ./bin/veqtl-mapper --het --bed data/phenotype.bed --job-number 1 --genes 10 --vcf data/genotype.vcf.gz --perm 10000,4

  const auto optsHet = new Opts(("./bin/veqtl-mapper --het --bed data/phenotype.bed --job-number 1 --genes 10 --vcf data/genotype.vcf.gz --perm 10000,4 --out " ~ testFile)
      .split);

  outFile = makeOut(optsHet);

  foreach (ref e; phenotype)
  {
    analyseData(e, permutations, outFile, optsHet, orderBuffer);
  }

  outFile.close;

  hash.start;
  put(hash, File(testFile).byChunk(1024));
  assert(toHexString(hash.finish) == "600B468E56A76D97B2B6CB1775A431F47ED4D65E");
  stderr.writeln("Passed: heterozygote test.");

  const auto optsNormal = new Opts(("./bin/veqtl-mapper --normal --bed data/phenotype.bed --job-number 1 --genes 10 --vcf data/genotype.vcf.gz --perm 10000,4 --out " ~ testFile)
      .split);

  outFile = makeOut(optsNormal);

  foreach (ref e; phenotype)
  {
    analyseData(e, permutations, outFile, optsNormal, orderBuffer);
  }

  outFile.close;

  hash.start;
  put(hash, File(testFile).byChunk(1024));
  assert(toHexString(hash.finish) == "194CF372ACA1414B3927238EEA35324E9C40BE16");
  stderr.writeln("Passed: normal test.");
}
