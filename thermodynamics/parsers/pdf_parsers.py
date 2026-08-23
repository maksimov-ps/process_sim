from __future__ import annotations
import shutil
import subprocess


class DockerNotAvailableError(RuntimeError):
    """Raised when Docker is not available on the system."""
    pass


class PDFParser(): 

    def __init__(self, 
                 verbose: bool = True): 

        self.verbose = verbose


    def _verify_docker_is_available(self) -> None: 
        if shutil.which("docker") is None:
            raise DockerNotAvailableError("Docker CLI is not installed or not available")

    
    def _check_if_container_exists(self, container_name: str | None) -> bool:

        if container_name is None:
            return False

        inspect_process = subprocess.run(
            ["docker", "container", "inspect", container_name],
            check=False,
            capture_output=True,
            text=True,
        )

        if self.verbose:
            if inspect_process.returncode == 0:
                print(f"Container '{container_name}' exists.")
            else:
                print(f"Could not find container '{container_name}'. Error: {inspect_process.stderr}")

        return inspect_process.returncode == 0


    def _start_existing_container(self, 
                                  container_name: str, 
                                  detach: bool = True,
                                  dry_run: bool = False) -> subprocess.CompletedProcess[str]:

        docker_command = ["docker", "start"]
        if not detach:
            docker_command.append("-a")
        docker_command.append(container_name)

        if self.verbose and not dry_run:
            print(f"Starting existing container '{container_name}'...")

        if dry_run:
            print(f"Dry run: {' '.join(docker_command)}")
            return docker_command

        return subprocess.run(
            docker_command,
            check=True,
            capture_output=True,
            text=True,
        )


    def _run_new_container(self,
                           image: str,
                           container_name: str | None = None,
                           detach: bool = True,
                           remove: bool = True,
                           dry_run: bool = False) -> subprocess.CompletedProcess[str]:

        docker_command: list[str] = ["docker", "run"]

        if remove:
            docker_command.append("--rm")
        if detach:
            docker_command.append("-d")
        if container_name:
            docker_command.extend(["--name", container_name])

        docker_command.append(image)

        if self.verbose and not dry_run:
            print(f"Running new container from image '{image}'...")

        if dry_run:
            print(f"Dry run: {' '.join(docker_command)}")
            return docker_command

        return subprocess.run(
            docker_command,
            check=True,
            capture_output=True,
            text=True,
        )


    def _initiate_docker_container(
        self,
        image: str,
        *,
        container_name: str | None = None,
        detach: bool = True,
        remove: bool = True,
        dry_run: bool = True,
    ) -> list[str] | subprocess.CompletedProcess[str]:

        """
        Initiates or reuses a Docker container.

        Parameters:
            image (str): The Docker image to run.
            container_name (str | None): The name of the container. If a container with this
                name already exists, this method starts it instead of running a new image.
            detach (bool): Whether to run the container in detached mode.
            remove (bool): Whether to remove the container after it exits when creating a new
                container with docker run.
            dry_run (bool): If True, return the Docker command instead of executing it.

        Returns:
            list[str] | subprocess.CompletedProcess[str]: The Docker command if dry_run is True,
            otherwise the result of subprocess.run.
        """

        if self._check_if_container_exists(container_name): 
            subprocess = self._start_existing_container(
                container_name=container_name, 
                detach=detach, 
                dry_run=dry_run
            )
            return subprocess


        subprocess = self._run_new_container(
            image=image,
            container_name=container_name,
            detach=detach,
            remove=remove,
            dry_run=dry_run,
        )
        return subprocess


class VLEDataPDFParser(PDFParser): 

    def __init__(self): 
        super().__init__()

        self._verify_docker_is_available()

        grobid = self._initiate_docker_container(
            image="grobid/grobid:0.9.0-full",
            container_name="pavels_grobid",
            detach=True,
            remove=True,
            dry_run=True,
        )

        pass 











