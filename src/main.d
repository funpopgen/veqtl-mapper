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

import core.stdc.stdlib : exit;
import std.conv : to;
import std.process : executeShell;
import std.range : enumerate;
import std.stdio : File, stderr, writeln;

import arg_parse : Opts;
import read_data : Phenotype, readBed, makeOut;
import run_analysis : analyseData;
import calculation : genPerms;

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
  pragma(msg, "VEQM");

  auto checkTabix = executeShell("command -v tabix");

  if (checkTabix.status != 0)
  {
    stderr.writeln("Error: tabix is not installed.");
    exit(1);
  }

  auto opts = new Opts(args.to!(string[]));

  if (opts.verbose)
  {
    stderr.writeln("Looking for ", opts.het ? "parent of origin" : "variance", " effects.");
  }

  auto phenotype = readBed(opts);

  if (opts.verbose)
  {
    stderr.writeln("Read ", phenotype.length, " phenotypes.");
  }

  auto permutations = genPerms(opts, phenotype[0].values.length);

  auto outFile = makeOut(opts);

  auto orderBuffer = new size_t[](phenotype[0].values.length);

  foreach (ref e; phenotype.enumerate)
  {
    if (opts.verbose)
    {
      stderr.writeln("Analysing gene ", e[1].geneName, " (", e[0] + 1,
          " out of ", phenotype.length, ").");
    }
    analyseData(e[1], permutations, outFile, opts, orderBuffer);
  }
}

unittest
{
  import core.stdc.stdlib : exit;
  import std.digest.sha;
  import std.file : exists, remove;
  import std.range : put;
  import std.stdio : stderr;

  if ("testtemp".exists)
  {
    stderr.writeln("Running unittests would overwrite testtemp file");
    exit(0);
  }

  scope (exit)
  {
    if ("testtemp".exists)
      "testtemp".remove;
  }

  // ./bin/VEQM --bed data/phenotype.bed --job-number 1 --genes 10 --vcf data/genotype_veqm.vcf.gz --perm 10000,4
  auto args = [
    "prog", "--bed", "data/phenotype.bed", "--job-number", "1", "--genes", "10",
    "--vcf", "data/genotype_veqm.vcf.gz", "--perm", "10000,4", "--out", "testtemp"
  ];

  auto opts = new Opts(args.to!(string[]));

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
  put(hash, File("testtemp").byChunk(1024));
  assert(toHexString(hash.finish) == "26F6F2F43BB3543FE7A8E060F0F9A088A74EACA1");
  stderr.writeln("Passed: variance test.");

  // ./bin/VEQM --het --bed data/phenotype.bed --job-number 1 --genes 10 --vcf data/genotype_veqm.vcf.gz --perm 10000,4

  opts.het = true;

  outFile = makeOut(opts);

  foreach (ref e; phenotype)
  {
    analyseData(e, permutations, outFile, opts, orderBuffer);
  }

  outFile.close;

  hash.start;
  put(hash, File("testtemp").byChunk(1024));
  assert(toHexString(hash.finish) == "600B468E56A76D97B2B6CB1775A431F47ED4D65E");
  stderr.writeln("Passed: heterozygote test.");

}
