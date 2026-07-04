#include <stdlib.h>
#include <string.h>
#include <sys/wait.h>
#include <unistd.h>

#include "minunit.h"

#include "phyc/mjson.h"

static const json_field schema[] = {
    {"sitepattern", JSON_REQUIRED, JSON_STRING},
    {"epsilon", JSON_OPTIONAL, JSON_NUMBER},
    {"invariant", JSON_OPTIONAL, JSON_BOOL},
    {"model", JSON_OPTIONAL, JSON_OBJECT | JSON_STRING},
    {"substitutionmodel", JSON_FORBIDDEN, JSON_ANY, "renamed to 'model'"},
};
static const size_t schema_n = sizeof(schema) / sizeof(schema[0]);

// Run json_validate in a forked child so its exit() does not kill the test
// runner. Returns the child's exit code (12 on a validation failure, 0 on
// success).
static int validate_status(const char* json) {
    // Flush first so the child does not inherit (and later re-flush on exit())
    // our buffered stdout, which would duplicate the harness' output.
    fflush(stdout);
    fflush(stderr);
    pid_t pid = fork();
    if (pid == 0) {
        // Silence diagnostics and harness noise so test output stays clean.
        freopen("/dev/null", "w", stderr);
        freopen("/dev/null", "w", stdout);
        json_node* root = create_json_tree(json);
        json_validate(root, schema, schema_n);
        json_free_tree(root);
        _exit(0);
    }
    int status = 0;
    waitpid(pid, &status, 0);
    return WIFEXITED(status) ? WEXITSTATUS(status) : -1;
}

char* test_valid() {
    const char* json =
        "{\"id\":\"s\",\"type\":\"sitemodel\",\"sitepattern\":\"sp\","
        "\"epsilon\":0.1,\"_comment\":\"ignored\"}";
    mu_assert(validate_status(json) == 0, "valid node should pass");
    return NULL;
}

char* test_missing_required() {
    const char* json = "{\"id\":\"s\",\"type\":\"sitemodel\",\"epsilon\":0.1}";
    mu_assert(validate_status(json) == 12, "missing required key should die");
    return NULL;
}

char* test_missing_id_type() {
    const char* json = "{\"sitepattern\":\"sp\"}";
    mu_assert(validate_status(json) == 12, "missing id/type should die");
    return NULL;
}

char* test_unknown_key() {
    const char* json =
        "{\"id\":\"s\",\"type\":\"sitemodel\",\"sitepattern\":\"sp\","
        "\"epsilonn\":0.1}";
    mu_assert(validate_status(json) == 12, "unknown key should die");
    return NULL;
}

char* test_forbidden_key() {
    const char* json =
        "{\"id\":\"s\",\"type\":\"sitemodel\",\"sitepattern\":\"sp\","
        "\"substitutionmodel\":\"m\"}";
    mu_assert(validate_status(json) == 12, "forbidden key should die");
    return NULL;
}

char* test_wrong_type() {
    // sitepattern declared JSON_STRING but given a number
    const char* json =
        "{\"id\":\"s\",\"type\":\"sitemodel\",\"sitepattern\":3}";
    mu_assert(validate_status(json) == 12, "wrong value type should die");
    return NULL;
}

char* test_object_or_string() {
    const char* ref =
        "{\"id\":\"s\",\"type\":\"sitemodel\",\"sitepattern\":\"sp\","
        "\"model\":\"m\"}";
    const char* obj =
        "{\"id\":\"s\",\"type\":\"sitemodel\",\"sitepattern\":\"sp\","
        "\"model\":{\"id\":\"m\",\"type\":\"x\"}}";
    mu_assert(validate_status(ref) == 0, "string reference should pass");
    mu_assert(validate_status(obj) == 0, "inline object should pass");
    return NULL;
}

// Run json_validate_xor on "model"/"function" in a forked child.
static int xor_status(const char* json) {
    fflush(stdout);
    fflush(stderr);
    pid_t pid = fork();
    if (pid == 0) {
        freopen("/dev/null", "w", stderr);
        freopen("/dev/null", "w", stdout);
        json_node* root = create_json_tree(json);
        json_validate_xor(root, "model", "function", NULL);
        json_free_tree(root);
        _exit(0);
    }
    int status = 0;
    waitpid(pid, &status, 0);
    return WIFEXITED(status) ? WEXITSTATUS(status) : -1;
}

char* test_xor() {
    const char* one = "{\"id\":\"s\",\"type\":\"t\",\"model\":\"m\"}";
    const char* none = "{\"id\":\"s\",\"type\":\"t\"}";
    const char* both =
        "{\"id\":\"s\",\"type\":\"t\",\"model\":\"m\",\"function\":\"f\"}";
    mu_assert(xor_status(one) == 0, "exactly one should pass");
    mu_assert(xor_status(none) == 12, "neither defined should die");
    mu_assert(xor_status(both) == 12, "both defined should die");
    return NULL;
}

char* test_required_getter() {
    const char* json = "{\"id\":\"s\",\"type\":\"t\",\"rate\":0.5}";
    json_node* root = create_json_tree(json);
    mu_assert(get_json_node_value_double_required(root, "rate") == 0.5,
              "required double getter wrong value");
    mu_assert(strcmp(get_json_node_value_string_required(root, "id"), "s") == 0,
              "required string getter wrong value");
    json_free_tree(root);
    return NULL;
}

char* all_tests() {
    char* message = NULL;
    mu_run_test(test_valid);
    mu_run_test(test_missing_required);
    mu_run_test(test_missing_id_type);
    mu_run_test(test_unknown_key);
    mu_run_test(test_forbidden_key);
    mu_run_test(test_wrong_type);
    mu_run_test(test_object_or_string);
    mu_run_test(test_xor);
    mu_run_test(test_required_getter);
    return NULL;
}

RUN_TESTS(all_tests)
