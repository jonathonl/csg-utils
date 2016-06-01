package CSG::Storage::SlotCtl::Command::list_pools;

use CSG::Storage::SlotCtl -command;
use CSG::Storage::Slots::DB;
use CSG::Base;

my $schema = CSG::Storage::Slots::DB->new();

sub opt_spec {
  return (
    ['project|p=s', 'List pools in a specific project [default: topmed]', {default => 'topmed'}],
  );
}

sub validate_args {
  my ($self, $opts, $args) = @_;

  if ($opts->{project}) {
    unless ($schema->resultset('Project')->find({name => $opts->{project}})) {
      $self->usage_error('Unable to locate project');
    }
  }
}

sub execute {
  my ($self, $opts, $args) = @_;

  my $pools = $schema->resultset('Pool');
  if ($opts->{project}) {
    $pools = $pools->search({'project.name' => $opts->{project}}, {join => 'project'});
  }

  for my $pool ($pools->all()) {
    say 'Name: ' . $pool->name;
    say "\tID: " . $pool->id;
    say "\tProject: " . $pool->project->name;
    say "\tHostname: " . $pool->hostname;
    say "\tPath: " . $pool->path;
    say "\tSpace Total: " . $pool->size_total;
    say "\tSpace Used: " . $pool->size_used;
    say '';
  }
}

1;

__END__

=head1

CSG::Storage::SlotCtl::Command::list_pools - List all defined pools
