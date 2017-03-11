#include <stdio.h>
#include <argp.h>

//Argument parsing
const char *argp_program_version = "mpcgs 0.2.0a";
static char doc[] = "mpcgs -- a parallel coalescent genealogy sampler.";

static struct argp_option options[] = {
    {"verbose", 'v', "LEVEL",  OPTION_ARG_OPTIONAL, "Set verbosity LEVEL" },
    {"quiet",   'q', 0, 0, "Produce minimal output" },
    {"gpu",     'g', 0, 0, "use GPU for theta estimation"},
    {"numiter", 'n', "INTEGER", 0, "run INTEGER estimation iterations"},
    {"chainlen", 'c', "INTEGER", 0, "generated INTEGER samples per iteration"},
    {"burninlen", 'b', "INTEGER", 0, "length of burn-in phase"},
    {"theta",   't', "VALUE", 0, "Set initial driving theta to VALUE" },
    {"input",   'i', "FILE", 0, "Sequence data contained in FILE (phylib format)"},
    { 0 }
};

static error_t parse_opt(int key, char *arg, struct argp_state *state)
{
  
    return(0);

}

static struct argp argp = {options, parse_opt, 0, doc};


///end argument parsing

int main(int argc, char *argv[])
{

    argp_parse(&argp, argc, argv, 0, 0, 0);

    printf("Theta estimate: 1.521\n");

    return(0);

}
