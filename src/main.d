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

  auto phenotype = readBed(opts);

  auto permutations = genPerms(opts, phenotype[0].values.length);

  auto outFile = makeOut(opts);

  auto orderBuffer = new size_t[](phenotype[0].values.length);

  foreach (ref e; phenotype)
  {
    analyseData(e, permutations, outFile, opts, orderBuffer);
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
  assert(toHexString(hash.finish) == "166F762354D571A6E44DDC7C04318C5FE9B11866");
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
  assert(toHexString(hash.finish) == "5317CB6C7EDE41C8B06986FF039B5DBA369859F3");
  stderr.writeln("Passed: heterozygote test.");

}
