Automatic n-dimensional recurrences over finite p-groups
with negative predecessor coordinates
Bogdan Dumitru, Mihai Prunescu

This folder contains the code, input images and certificates
needed to reproduce its reported computations. All paths below are relative
to this folder. The programs require no files from the surrounding workspace.

QUICK START

Run all computational checks:

    python3 reproduce.py

Or, with make installed:

    make check

The programs use Python 3 and its standard library, with no pip dependencies.
They were tested with Python 3.14.7. Run with assertions enabled (no -O or
-OO). A successful run ends with:

    All reported computational checks passed.

CODE AND THE RESULTS IT CHECKS

1. Binary automata and block substitutions

   Entry point: research/certify_original_examples.py
   Shared routines: research/automatic_recurrences.py, research/png_read.py
   and examples/compute_examples.py.

   Checks both binary examples in "Two certified binary examples". It verifies
   the recurrence and boundary conditions at every reachable product state,
   minimizes the automata, checks numerical and colored substitutions, and
   compares their outputs with the four supplied PNG images in original/.

   Expected minimal automaton sizes (least digits first / most digits first):
       negative predecessor example: 53 / 107
       five nonnegative predecessors: 254 / 710
   Expected reachable verification states: 686 and 552493.
   Expected numerical block counts: 44 and 204.
   Expected colored substitution label counts: 44 and 212.

   The main runner also reconstructs the two machines from finite samples and
   recovers the rejected candidates' witnesses (2,2), (4,4), and (52,16).

2. Nonabelian examples and figures

   Entry point: examples/nonabelian_examples.py

   Recomputes both group recurrences, verifies the group maps, compares the
   two multiplication orders, and regenerates the figures and counts in
   "Nonabelian examples" and "Exact computations on the displayed squares".

   On the 243 by 243 square, expected counts are:
       39 distinct values in the group of order 81;
       all 27 values in the group of order 27;
       17737 positions changed by reversing multiplication.

   The runner compares the regenerated pixels with the four paper figures.

3. Delannoy counts and automata

   Entry point: research/verify_mod27.py
   Shared routines: examples/compute_examples.py.

   Checks the 1362-state automaton modulo 27, all 12258 polynomial transition
   identities, reachability and minimality. It verifies the integer-vector
   identity M^4 (M-8I)^3 v = 0 and the initial counts 6165, 54384, 477888.
   Together these certify the paper's formula for every r >= 4:

       N_27(r) = 8^(r-4) (18r^2 + 471r + 3993).

   It also compares the automaton with 59049 residues evaluated directly from
   the Delannoy recurrence. The main runner constructs the machines modulo
   3 and 9 and checks their minimal state counts, 3 and 58.

Run an individual entry point with python3 followed by its path. The other
Python files are shared modules. All three entry points are run by reproduce.py.

DATA AND OUTPUTS

research/generated/
    Supplied transition tables, polynomial states, block substitutions,
    discovery records and certificates. The programs regenerate the reports
    and certificates. reproduce.py writes reproduction_results.json here.

examples/generated/nonabelian/
    Regenerated images, small_arrays.json and verification.json.

original/
    Four input images for the two binary examples, including their colored
    versions. These images are required by the binary verifier.

images/
    The four PNG files included in the paper. They are also the reference
    images checked by reproduce.py.

The certificates for automata and the modulo-27 formula verify identities
that apply at all indices. Counts for the nonabelian examples concern the
specified finite square. The general automaticity theorems are proved in
the paper.

