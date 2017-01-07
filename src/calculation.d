module calculation;

import arg_parse : Opts;
import std.algorithm : sum;
import std.conv : to;
import std.exception : enforce;
import std.math : fabs, sqrt;

version (unittest)
{
  import std.math : approxEqual;
}

class VarianceException : Exception
{
  //thrown if variable is constant
  pure nothrow this(string s)
  {
    super(s);
  }
}

pure nothrow extern (C)
{
  //call GSL to calculate P values from T statistics
  double gsl_cdf_tdist_P(double x, double nu);
}

extern (C)
{
  int mleBeta(double* pval, size_t nPerm, double* alpha, double* beta);
}

@system unittest
{
  // Checks GSL gives right p value for t statistic
  assert(approxEqual(gsl_cdf_tdist_P(-1.6, 7), 0.07681585));
}

pure ref T[] rank(T)(ref T[] rankArray, size_t[] orderIndex)
{
  //ranks array, giving ties mean rank
  import std.algorithm : makeIndex;

  immutable size_t len = rankArray.length;
  makeIndex!("a < b")(rankArray, orderIndex);

  T sumrank = 0.0;
  size_t dupcount = 0;
  T avgrank;

  foreach (i, ref e; orderIndex)
  {
    sumrank += i;
    dupcount++;
    if (i == (len - 1) || rankArray[e] != rankArray[orderIndex[i + 1]])
    {
      avgrank = sumrank / dupcount + 1;
      foreach (ref j; orderIndex[(i - dupcount + 1) .. (i + 1)])
        rankArray[j] = avgrank;
      sumrank = 0;
      dupcount = 0;
    }
  }
  return rankArray;
}

@safe unittest
{
  //Simple test of ranking with ties
  auto vector = [10.0, 9, 2, 9, 3];
  auto orderBuffer = new size_t[](5);

  assert(rank(vector, orderBuffer) == [5, 3.5, 1, 3.5, 2]);
}

pure ref T[] rank_discrete(T)(ref T[] rankArray)
{
  import std.conv : to;

  //easier to rank discrete genotypes
  T[3] countGenotypes = 0;
  T[3] rankGenotypes;

  //count numbers of 0, 1, 2 alleles
  foreach (ref e; rankArray)
    countGenotypes[e.to!size_t]++;
  rankGenotypes[0] = (countGenotypes[0] + 1) / 2;
  rankGenotypes[1] = (2 * countGenotypes[0] + countGenotypes[1] + 1) / 2;
  rankGenotypes[2] = (2 * countGenotypes[0] + 2 * countGenotypes[1] + countGenotypes[2] + 1) / 2;

  foreach (ref e; rankArray)
    e = rankGenotypes[e.to!size_t];
  return rankArray;
}

@safe unittest
{
  //Ranking of discrete genotypes
  import std.algorithm : map;
  import std.array : array;
  import std.random : uniform;
  import std.range : iota;

  auto vector = [0.0, 1, 0, 2, 2, 2, 1];
  assert(rank_discrete(vector) == [1.5, 3.5, 1.5, 6, 6, 6, 3.5]);

  auto ranVector = iota(20).map!(a => uniform(0, 3).to!double).array;
  auto ranVector2 = ranVector.dup;
  auto orderBuffer = new size_t[](20);
  assert(rank_discrete(ranVector) == rank(ranVector2, orderBuffer));
}

pure void transform(ref double[] vector)
{
  //transforms array so mean =0 sum of squares = 1
  int n = 0;
  double mean = 0;
  double M2 = 0;
  double delta;

  foreach (ref e; vector)
  {
    n++;
    delta = e - mean;
    mean += delta / n;
    M2 += delta * (e - mean);
  }

  enforce(M2 != 0, new VarianceException(""));

  M2 = sqrt(M2);

  foreach (ref e; vector)
    e = (e - mean) / M2;
}

@system unittest
{
  //Checks that transform works on randomly generated vector
  import std.algorithm : reduce;
  import std.random : uniform;

  auto x = new double[](10);
  foreach (ref e; x)
    e = uniform(0.0, 10.0);

  transform(x);
  auto mean = 0.0.reduce!((a, b) => a + b)(x);

  assert(approxEqual(mean, 0.0));
  assert(approxEqual(0.0.reduce!((a, b) => a + (b - mean) * (b - mean))(x), 1));
}

pure nothrow double[2] correlation(ref double[] vector1, ref double[] vector2)
{
  //calculates correlation, t stat and p value for two arrays
  import std.numeric : dotProduct;

  double[2] results;
  results[0] = dotProduct(vector1, vector2);
  results[1] = gsl_cdf_tdist_P(
      -fabs(results[0] * sqrt((vector1.length - 2) / (1 - results[0] * results[0]))),
      vector1.length - 2) * 2;
  return results;
}

double corPvalue(double cor, size_t nInd)
{
  return gsl_cdf_tdist_P(-fabs(cor * sqrt((nInd - 2) / (1 - cor * cor))), nInd - 2) * 2;
}

@system unittest
{
  //Check correlation of phenotype with 3rd row genotype against estimates from R
  auto corFromR = [-0.2863051, 0.4225695];

  auto genotype = [0.115, 2, 0.0964, 1, 1, 1, 0, 1, 0, 0.0563];
  auto phen = [
    -1.3853088072, -0.785797093643, 1.14540423638, -0.785797093643, 1.03820492508,
    -1.25652676836, -0.787662180447, -2.05355237841, -0.245457234103, 1.14277217712
  ];

  auto orderBuffer = new size_t[](phen.length);

  transform(rank(phen, orderBuffer));
  transform(rank(genotype, orderBuffer));
  auto cor = correlation(genotype, phen);

  assert(approxEqual(cor[0], corFromR[0]));
  assert(approxEqual(cor[1], corFromR[1]));
}

size_t[] genPerms(const Opts opts, size_t nInd)
{
  import std.array : array;
  import std.random : randomShuffle, rndGen;
  import std.range : chunks, cycle, iota, take;

  if (opts.perms.length > 1)
    rndGen.seed(opts.perms[1]);

  size_t[] outPerm = iota(nInd).cycle.take(opts.perms[0] * nInd).array;

  foreach (ref perm; chunks(outPerm, nInd))
    randomShuffle(perm);

  return outPerm;

}

double[2] betaParameters(ref double[] minPval)
{
  immutable auto mean = minPval.sum / minPval.length;
  auto variance = 0.0;
  foreach (ref e; minPval)
    variance += (e - mean) * (e - mean);
  variance /= minPval.length;
  //Estimate shape1 & shape2
  auto alpha = mean * (mean * (1 - mean) / variance - 1);
  auto beta = alpha * (1 / mean - 1);

  immutable auto alphaCopy = alpha;
  immutable auto betaCopy = beta;

  immutable auto success = mleBeta(minPval.ptr, minPval.length, &alpha, &beta);

  if (success == 0)
  {
    return [alphaCopy, betaCopy];
  }

  return [alpha, beta];
}
