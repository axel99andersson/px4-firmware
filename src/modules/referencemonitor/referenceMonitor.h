/****************************************************************************
 *
 *   Copyright (c) 2018 PX4 Development Team. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 * 3. Neither the name PX4 nor the names of its contributors may be
 *    used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
 * OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 * AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 * ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 ****************************************************************************/

 #pragma once

 #include <px4_platform_common/module.h>
 #include <px4_platform_common/module_params.h>

 extern "C" __EXPORT int referenceMonitor_main(int argc, char *argv[]);
 class ReferenceMonitor : public ModuleBase<ReferenceMonitor>, public ModuleParams
 {
 public:
	 ReferenceMonitor(int example_param, bool example_flag);

	 virtual ~ReferenceMonitor() = default;

	 /** @see ModuleBase */
	 static int task_spawn(int argc, char *argv[]);

	 /** @see ModuleBase */
	 static ReferenceMonitor *instantiate(int argc, char *argv[]);

	 /** @see ModuleBase */
	 static int custom_command(int argc, char *argv[]);

	 /** @see ModuleBase */
	 static int print_usage(const char *reason = nullptr);

	 /** @see ModuleBase::run() */
	 void run() override;

	 /** @see ModuleBase::print_status() */
	 int print_status() override;

 private:

	 /**
	  * Check for parameter changes and update them if needed.
	  * @param parameter_update_sub uorb subscription to parameter_update
	  * @param force for a parameter update
	  */
	 void parameters_update(int parameter_update_sub, bool force = false);


	 DEFINE_PARAMETERS(
		 (ParamInt<px4::params::SYS_AUTOSTART>) _sys_autostart,   /**< example parameter */
		 (ParamInt<px4::params::SYS_AUTOCONFIG>) _sys_autoconfig  /**< another parameter */
	 )
 };
